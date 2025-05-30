use std::pin::Pin;

use bytes::Bytes;
use futures::{Stream, TryStreamExt, stream};

use super::{Client, Error, response::ticket::BlockUrl};

pub(crate) fn chunks<'a>(
    client: &'a Client,
    urls: &'a [BlockUrl],
) -> impl Stream<Item = crate::Result<Bytes>> + 'a {
    Box::pin(
        stream::try_unfold((client, urls, 0), |(client, urls, i)| async move {
            match urls.get(i) {
                Some(url) => {
                    let st = resolve_data(client, url).await;
                    Ok(Some((st, (client, urls, i + 1))))
                }
                None => Ok(None),
            }
        })
        .try_flatten(),
    )
}

async fn resolve_data(
    client: &Client,
    block_url: &BlockUrl,
) -> Pin<Box<dyn Stream<Item = crate::Result<Bytes>>>> {
    use base64::prelude::{BASE64_STANDARD, Engine as _};

    const DELIMITER: &str = ";base64,";

    let url = block_url.url();

    if url.scheme() == "data" {
        // _Htsget retrieval API spec v1.3.0_ § "Inline data block URIs": "client should ignore the
        // media type (if any), treating the payload as a partial blob."
        if let Some((_, encoded_data)) = url.as_str().split_once(DELIMITER) {
            match BASE64_STANDARD.decode(encoded_data) {
                Ok(data) => Box::pin(stream::once(async { Ok(Bytes::from(data)) })),
                Err(e) => Box::pin(stream::once(async { Err(Error::Decode(e)) })),
            }
        } else {
            Box::pin(stream::once(async { Err(Error::InvalidDataUrl) }))
        }
    } else {
        let mut request = client.http_client().get(url.clone());

        for (key, value) in block_url.headers() {
            request = request.header(key, value);
        }

        match request.send().await.map_err(Error::Request) {
            Ok(response) => Box::pin(response.bytes_stream().map_err(Error::Request)),
            Err(e) => Box::pin(stream::once(async { Err(e) })),
        }
    }
}
