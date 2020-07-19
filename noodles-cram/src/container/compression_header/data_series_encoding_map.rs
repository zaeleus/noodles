pub mod data_series;

pub use self::data_series::DataSeries;

use super::Encoding;

/// A container compression header data series encoding map.
#[derive(Debug, Default)]
pub struct DataSeriesEncodingMap {
    map: [Option<Encoding>; DataSeries::LEN],
}

impl DataSeriesEncodingMap {
    /// Inserts a data series-encoding pair into the map.
    pub fn insert(&mut self, data_series: DataSeries, encoding: Encoding) {
        let i = data_series as usize;
        self.map[i] = Some(encoding);
    }

    /// Returns a reference to the value corresponding to the given data series.
    pub fn get(&self, data_series: &DataSeries) -> Option<&Encoding> {
        let i = *data_series as usize;
        self.map[i].as_ref()
    }
}
