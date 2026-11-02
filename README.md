# huffman_zandronum

A small, dependency-free Rust implementation of the Huffman codec used by Zandronum.

> **Compatibility:** This crate implements a fixed Zandronum-compatible Huffman tree. It is not a general-purpose Huffman compressor: encoded data is intended for software that uses the same protocol format and tree.

## Features

- Encode and decode byte slices using the built-in Zandronum-compatible tree.
- Reuse caller-owned `Vec<u8>` buffers to avoid allocating a new output vector for each call.
- Set a maximum output size for each operation.
- Match reference C++ codec vectors, including raw fallback, padding, and binary data.
- No third-party dependencies.

## Installation

Add it to your project:

```toml
[dependencies]
huffman_zandronum = "0.1.0"
```

Or run:

```sh
cargo add huffman_zandronum
```

## Quick start

```rust
use huffman_zandronum::{decode, encode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let message = b"Zandronum protocol data";
    let mut encoded = Vec::new();
    let mut decoded = Vec::new();

    // `max_len` limits the encoded bytes, including the framing byte.
    let encoded_len = encode(message, &mut encoded, 1024)?;

    // Decode into a separate reusable buffer. The limit applies to decoded bytes.
    let decoded_len = decode(&encoded[..encoded_len], &mut decoded, 1024)?;

    assert_eq!(&decoded[..decoded_len], message);
    Ok(())
}
```

Both functions clear the supplied output vector before writing. The returned `usize` is the number of valid bytes in that vector; use the returned length when passing a buffer onward. The functions do not require the output vector to be empty, so the same vector can be reused across calls.

## API

### `encode(input, output, max_len)`

Encodes `input` into `output` and returns the encoded length. `max_len` limits the encoded bytes, including the padding byte. If the Huffman candidate exceeds this limit while being built, the function returns `EncodeError::OutputTooLarge`.

The codec selects raw framing when the Huffman representation would be larger than the original data plus the raw marker. Consequently, output is not guaranteed to be smaller than input.

### `decode(input, output, max_len)`

Decodes a framed byte slice into `output` and returns the decoded length. `max_len` limits the decoded data, not the encoded input size. If the decoded result would exceed the limit, the function returns `DecodeError::OutputTooLarge`.

An empty input decodes to an empty output. A one-byte raw marker (`0xFF`) likewise represents an empty raw payload.

### Errors and buffer state

The only current error variant for either operation is `OutputTooLarge`. The output buffer is cleared at the start of every call. On an error, it may contain a partial result; treat its contents as invalid unless the call returned `Ok`.

## Format details

The wire framing and padding rules are described in [docs/format.md](docs/format.md). This is a fixed protocol codec and does not expose a way to select a tree or configure code lengths.

## Development

Run the test suite with:

```sh
cargo test
```

The tests include golden byte vectors from the original C++ codec as well as encode/decode round trips.
