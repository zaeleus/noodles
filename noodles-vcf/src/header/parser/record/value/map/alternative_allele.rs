use std::{error, fmt};

use crate::{
    header::record::value::{
        map::{
            self,
            alternative_allele::{tag, Tag},
            AlternativeAllele, OtherFields,
        },
        Map,
    },
    record::alternate_bases::allele::{symbol, Symbol},
};

#[derive(Clone, Debug, Eq, PartialEq)]
enum ParseErrorKind {
    InvalidMap(super::ParseError),
    InvalidField(super::field::ParseError),
    MissingId,
    InvalidId(symbol::ParseError),
    MissingDescription,
    DuplicateTag(Tag),
}

/// An error returned when a VCF header record alternative allele map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError {
    id: Option<Symbol>,
    kind: ParseErrorKind,
}

impl ParseError {
    fn new(id: Option<Symbol>, kind: ParseErrorKind) -> Self {
        Self { id, kind }
    }

    pub(crate) fn id(&self) -> Option<&Symbol> {
        self.id.as_ref()
    }
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match &self.kind {
            ParseErrorKind::InvalidMap(e) => Some(e),
            ParseErrorKind::InvalidField(e) => Some(e),
            ParseErrorKind::InvalidId(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ParseErrorKind::InvalidMap(_) => write!(f, "invalid map"),
            ParseErrorKind::InvalidField(_) => write!(f, "invalid field"),
            ParseErrorKind::MissingId => write!(f, "missing ID"),
            ParseErrorKind::InvalidId(_) => write!(f, "invalid ID"),
            ParseErrorKind::MissingDescription => write!(f, "missing description"),
            ParseErrorKind::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub fn parse_alternative_allele(
    src: &mut &[u8],
) -> Result<(Symbol, Map<AlternativeAllele>), ParseError> {
    super::consume_prefix(src).map_err(|e| ParseError::new(None, ParseErrorKind::InvalidMap(e)))?;

    let mut id = None;
    let mut description = None;
    let mut other_fields = OtherFields::new();

    while let Some((raw_key, raw_value)) = super::split_field(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidField(e)))?
    {
        match Tag::from(raw_key) {
            tag::ID => {
                parse_id(&raw_value, &id).and_then(|v| try_replace(&mut id, &None, tag::ID, v))?;
            }
            tag::DESCRIPTION => {
                try_replace(&mut description, &id, tag::DESCRIPTION, raw_value.into())?;
            }
            Tag::Other(t) => try_insert(&mut other_fields, &id, t, raw_value.into())?,
        }
    }

    super::consume_suffix(src)
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidMap(e)))?;

    let id = id.ok_or_else(|| ParseError::new(None, ParseErrorKind::MissingId))?;
    let description = description
        .ok_or_else(|| ParseError::new(Some(id.clone()), ParseErrorKind::MissingDescription))?;

    Ok((
        id,
        Map {
            inner: AlternativeAllele { description },
            other_fields,
        },
    ))
}

fn parse_id(s: &str, id: &Option<Symbol>) -> Result<Symbol, ParseError> {
    s.parse()
        .map_err(|e| ParseError::new(id.clone(), ParseErrorKind::InvalidId(e)))
}

fn try_replace<T>(
    option: &mut Option<T>,
    id: &Option<Symbol>,
    tag: Tag,
    value: T,
) -> Result<(), ParseError> {
    if option.replace(value).is_none() {
        Ok(())
    } else {
        Err(ParseError::new(
            id.clone(),
            ParseErrorKind::DuplicateTag(tag),
        ))
    }
}

fn try_insert(
    other_fields: &mut OtherFields<tag::Standard>,
    id: &Option<Symbol>,
    tag: map::tag::Other<tag::Standard>,
    value: String,
) -> Result<(), ParseError> {
    use indexmap::map::Entry;

    match other_fields.entry(tag) {
        Entry::Vacant(entry) => {
            entry.insert(value);
            Ok(())
        }
        Entry::Occupied(entry) => {
            let (t, _) = entry.remove_entry();
            Err(ParseError::new(
                id.clone(),
                ParseErrorKind::DuplicateTag(Tag::Other(t)),
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_filter() {
        use symbol::{structural_variant::Type, StructuralVariant};

        let mut src = &br#"<ID=DEL,Description="Deletion">"#[..];

        let id = Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion));
        let map = Map::<AlternativeAllele>::new("Deletion");
        let expected = (id, map);

        assert_eq!(parse_alternative_allele(&mut src), Ok(expected));
    }
}
