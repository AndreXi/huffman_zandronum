# Encoded byte format

This crate implements a fixed Huffman tree and framing compatible with the Zandronum codec. The format is byte-oriented: each input byte is one symbol, so input and decoded output may contain any value from `0x00` through `0xFF`. It is not UTF-8-specific.

## Frame layout

Every encoded frame starts with one framing byte.

- **Huffman frame (`0x00`–`0x07`):** the first byte is the number of unused padding bits at the end of the final payload byte. The remaining bytes contain the Huffman bitstream.
- **Raw frame (`0xFF`):** the remaining bytes are copied and are not Huffman-coded. The encoder chooses this form when the Huffman form would be longer than the input plus the marker.

The encoder emits an empty input as `[0x00]`: zero padding bits and no payload. A raw empty payload, `[0xFF]`, also decodes to empty output.

Values `0x08`–`0xFE` are not emitted as framing bytes by this encoder. The decoder treats a non-`0xFF` first byte as a padding count; malformed padding that exceeds the available payload bits produces an empty result rather than a format error.

## Bit order and padding

Huffman codes are traversed from their most-significant code bit to their least-significant bit. The stream is packed into bytes and each payload byte is bit-reversed for the protocol's on-wire ordering. The decoder reverses that operation before traversing the tree.

If the final payload byte is only partially used, the unused low-order positions in the traversal stream are padding. The leading framing byte records how many bits to ignore. Padding is not decoded as symbols. When the stream ends exactly on a byte boundary, the padding count is zero.

## Raw fallback and size

The raw marker costs one byte, so raw framing uses `input.len() + 1` bytes. The encoder switches to raw framing only when the Huffman representation is strictly larger than that size; ties stay Huffman-coded. This behavior is part of byte-for-byte compatibility.

The `max_len` argument to `encode` caps the output buffer size, including the framing byte. Encoding constructs the Huffman candidate before deciding whether raw fallback is preferable, so it can return `OutputTooLarge` if that candidate crosses the limit—even when the eventual raw frame would fit. The `max_len` argument to `decode` caps the decoded payload size. Both functions can return `OutputTooLarge`; callers handling untrusted or large messages should set limits appropriate to their application.

## Compatibility scope

The Huffman tree is embedded in the crate. The implementation's golden tests compare its output against reference C++ codec vectors. Use the codec only where the peer expects this same tree, bit order, and frame format like a Zandronum server.
