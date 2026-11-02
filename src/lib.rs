//! Zandronum-compatible Huffman encoding and decoding for byte slices.
//!
//! The codec uses a fixed Huffman tree and the protocol's framing format; it is
//! not a configurable or general-purpose Huffman implementation. Encoded and
//! decoded data may contain arbitrary bytes, not just text.
//!
//! Reuse a caller-owned output buffer to avoid allocating a new vector for each
//! operation:
//!
//! ```
//! use huffman_zandronum::{decode, encode};
//!
//! let mut encoded = Vec::new();
//! let mut decoded = Vec::new();
//! let message = b"protocol data";
//!
//! let encoded_len = encode(message, &mut encoded, 1024)?;
//! let decoded_len = decode(&encoded[..encoded_len], &mut decoded, 1024)?;
//! assert_eq!(&decoded[..decoded_len], message);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! See the package README for installation and [`encode`] / [`decode`] for
//! operation details.

const COMPATIBLE_HUFFMAN_TREE: &[u8] = &[
    0, 0, 0, 1, 128, 0, 0, 0, 3, 38, 34, 2, 1, 80, 3, 110, 144, 67, 0, 2, 1, 74, 3, 243, 142, 37,
    2, 3, 124, 58, 182, 0, 0, 1, 36, 0, 3, 221, 131, 3, 245, 163, 1, 35, 3, 113, 85, 0, 1, 41, 1,
    77, 3, 199, 130, 0, 1, 206, 3, 185, 153, 3, 70, 118, 0, 3, 3, 5, 0, 0, 1, 24, 0, 2, 3, 198,
    190, 63, 2, 3, 139, 186, 75, 0, 1, 44, 2, 3, 240, 218, 56, 3, 40, 39, 0, 0, 2, 2, 3, 244, 247,
    81, 65, 0, 3, 9, 125, 3, 68, 60, 0, 0, 1, 25, 3, 191, 138, 3, 86, 17, 0, 1, 23, 3, 220, 178, 2,
    3, 165, 194, 14, 1, 0, 2, 2, 0, 0, 2, 1, 208, 3, 150, 157, 181, 1, 222, 2, 3, 216, 230, 211, 0,
    2, 2, 3, 252, 141, 10, 42, 0, 2, 3, 134, 135, 104, 1, 103, 3, 187, 225, 95, 32, 0, 0, 0, 0, 0,
    0, 1, 57, 1, 61, 3, 183, 237, 0, 0, 3, 233, 234, 3, 246, 203, 2, 3, 250, 147, 79, 1, 129, 0, 1,
    7, 3, 143, 136, 1, 20, 3, 179, 148, 0, 0, 0, 3, 28, 106, 3, 101, 87, 1, 66, 0, 3, 180, 219, 3,
    227, 241, 0, 1, 26, 1, 251, 3, 229, 214, 3, 54, 69, 0, 0, 0, 0, 0, 3, 231, 212, 3, 156, 176, 3,
    93, 83, 0, 3, 96, 253, 3, 30, 13, 0, 0, 2, 3, 175, 254, 94, 3, 159, 27, 2, 1, 8, 3, 204, 226,
    78, 0, 0, 0, 3, 107, 88, 1, 31, 3, 137, 169, 2, 2, 3, 215, 145, 6, 4, 1, 127, 0, 1, 99, 3, 209,
    217, 0, 3, 213, 238, 3, 177, 170, 1, 132, 0, 0, 0, 2, 3, 22, 12, 114, 2, 2, 3, 158, 197, 97,
    45, 0, 1, 46, 1, 112, 3, 174, 249, 0, 3, 224, 102, 2, 3, 171, 151, 193, 0, 0, 0, 3, 15, 16, 3,
    2, 168, 1, 49, 3, 91, 146, 0, 1, 48, 3, 173, 29, 0, 3, 19, 126, 3, 92, 242, 0, 0, 0, 0, 0, 0,
    3, 205, 192, 2, 3, 235, 149, 255, 2, 3, 223, 184, 248, 0, 0, 3, 108, 236, 3, 111, 90, 2, 3,
    117, 115, 71, 0, 0, 3, 11, 50, 0, 3, 188, 119, 1, 122, 3, 167, 162, 1, 160, 1, 133, 3, 123, 21,
    0, 0, 2, 1, 59, 2, 3, 155, 154, 98, 43, 0, 3, 76, 51, 2, 3, 201, 116, 72, 2, 0, 2, 3, 109, 100,
    121, 2, 3, 195, 232, 18, 1, 0, 2, 0, 1, 164, 2, 3, 120, 189, 73, 0, 1, 196, 3, 239, 210, 3, 64,
    62, 89, 0, 0, 1, 33, 2, 3, 228, 161, 55, 2, 3, 84, 152, 47, 0, 0, 2, 3, 207, 172, 140, 3, 82,
    166, 0, 3, 53, 105, 1, 52, 3, 202, 200,
];

const LEAF_COUNT: usize = 256;
const BRANCH_COUNT: usize = LEAF_COUNT - 1;
const NODE_COUNT: usize = LEAF_COUNT + BRANCH_COUNT;
const REVERSE_BITS: bool = true;
const EXPANDABLE: bool = false;

#[derive(Clone, Copy)]
enum Node {
    Branch { left: u16, right: u16 },
    Leaf(u8),
}

/// The precomputed nodes and symbol codes for the built-in Huffman tree.
struct CodecData {
    nodes: [Node; NODE_COUNT],
    codes: [(u16, u8); LEAF_COUNT],
}

const fn build_child(
    data: &[u8],
    descriptor_index: &mut usize,
    descriptor: u8,
    direction_bit: u8,
    code: u16,
    bit_count: u8,
    nodes: &mut [Node; NODE_COUNT],
    next_node: &mut usize,
    codes: &mut [(u16, u8); LEAF_COUNT],
) -> u16 {
    // Appends the direction bit "0" left, "1" right
    let child_code: u16 = (code << 1) | direction_bit as u16;
    // Increments the depth counter
    let child_bits = bit_count + 1;

    // Bitmasks the descriptor to know if we have a branch or a leaf
    if descriptor & (1u8 << direction_bit) == 0 {
        build_node(
            data,
            descriptor_index,
            child_bits,
            child_code,
            nodes,
            next_node,
            codes,
        )
    } else {
        // Is a leaf
        let index = *next_node;

        // Move to the next node
        *next_node += 1;

        // Get the value and move the descriptor
        let value = data[*descriptor_index];
        *descriptor_index += 1;

        // Write the leaf node
        nodes[index] = Node::Leaf(value);
        codes[value as usize] = (child_code, child_bits);

        // We can't use safer methods in const fn
        // But still it should not be greater than u16
        index as u16
    }
}

const fn build_node(
    data: &[u8],
    descriptor_index: &mut usize,
    bit_count: u8,
    code: u16,
    nodes: &mut [Node; NODE_COUNT],
    next_node: &mut usize,
    codes: &mut [(u16, u8); LEAF_COUNT],
) -> u16 {
    // Set the index of the node and increment the index
    let index = *next_node;
    *next_node += 1;

    // Also with the descriptor set and increment
    let descriptor = data[*descriptor_index];
    *descriptor_index += 1;

    let left = build_child(
        data,
        descriptor_index,
        descriptor,
        0,
        code,
        bit_count,
        nodes,
        next_node,
        codes,
    );

    let right = build_child(
        data,
        descriptor_index,
        descriptor,
        1,
        code,
        bit_count,
        nodes,
        next_node,
        codes,
    );

    // Write the branch
    nodes[index] = Node::Branch { left, right };

    index as u16
}

const fn build_codec() -> CodecData {
    let mut nodes = [Node::Leaf(0); NODE_COUNT];
    let mut codes = [(0u16, 0u8); LEAF_COUNT];
    let mut descriptor_index: usize = 0;
    let mut next_node: usize = 0;

    build_node(
        COMPATIBLE_HUFFMAN_TREE,
        &mut descriptor_index,
        0,
        0,
        &mut nodes,
        &mut next_node,
        &mut codes,
    );

    CodecData {
        nodes: nodes,
        codes: codes,
    }
}

static HUFFMAN_ZANDRONUM: CodecData = build_codec();

/// An error returned while decoding.
#[derive(Debug)]
pub enum DecodeError {
    /// The decoded output would exceed the `max_len` limit.
    OutputTooLarge,
}

impl std::fmt::Display for DecodeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OutputTooLarge => formatter.write_str("decoded output exceeds max_len"),
        }
    }
}

impl std::error::Error for DecodeError {}

/// Decodes a Zandronum-framed byte slice into a reusable output buffer.
///
/// Clears `output` before decoding and returns the number of decoded bytes. The
/// decoded result is limited to `max_len`; the input frame itself is not size
/// limited. If decoding exceeds the limit, this returns
/// [`DecodeError::OutputTooLarge`] and `output` may contain partial data.
///
/// An empty input and a raw frame containing only its `0xFF` marker both decode
/// to empty output. See the package's `docs/format.md` for framing details.
///
/// # Examples
///
/// ```
/// # use huffman_zandronum::decode;
/// let mut output = Vec::new();
/// let len = decode(&[0xFF, b'h', b'i'], &mut output, 2)?;
/// assert_eq!(&output[..len], b"hi");
/// # Ok::<(), huffman_zandronum::DecodeError>(())
/// ```
pub fn decode(input: &[u8], output: &mut Vec<u8>, max_len: usize) -> Result<usize, DecodeError> {
    output.clear();
    let Some((&padding, payload)) = input.split_first() else {
        return Ok(0);
    };

    // Checks for unencoded signals.
    if padding == 0xFF {
        if payload.len() <= max_len {
            output.extend_from_slice(payload);
            return Ok(payload.len());
        } else {
            return Err(DecodeError::OutputTooLarge);
        }
    }

    // How many bits we can consume, we know we have at least 1 element here.
    // If we can't consume it's all padding.
    let mut bits_available: usize = match ((input.len() - 1) << 3).checked_sub(padding.into()) {
        Some(bits) => bits,
        None => return Ok(0),
    };

    let mut node: Node = HUFFMAN_ZANDRONUM.nodes[0];

    for &item in payload {
        let byte = if REVERSE_BITS {
            item.reverse_bits()
        } else {
            item
        };

        for i in (0..8).rev() {
            if bits_available == 0 {
                return Ok(output.len());
            }

            let bit = byte >> i & 1;
            bits_available -= 1;

            if let Node::Branch { left, right } = node {
                let child_index: usize = if bit == 0 { left } else { right }.into();
                node = HUFFMAN_ZANDRONUM.nodes[child_index];
            }

            if let Node::Leaf(value) = node {
                if output.len() < max_len {
                    output.push(value);
                } else {
                    return Err(DecodeError::OutputTooLarge);
                }
                node = HUFFMAN_ZANDRONUM.nodes[0];
            };
        }
    }

    Ok(output.len())
}

/// An error returned while encoding.
#[derive(Debug)]
pub enum EncodeError {
    /// The encoded output would exceed the `max_len` limit.
    OutputTooLarge,
}

impl std::fmt::Display for EncodeError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::OutputTooLarge => formatter.write_str("encoded output exceeds max_len"),
        }
    }
}

impl std::error::Error for EncodeError {}

/// Encodes a byte slice using the built-in Zandronum-compatible Huffman tree.
///
/// Clears `output` before encoding and returns the number of encoded bytes. The
/// `max_len` limit includes the framing byte. If the limit is exceeded, this
/// returns [`EncodeError::OutputTooLarge`] and `output` may contain partial
/// data. The encoder uses raw framing when the completed Huffman representation
/// is larger than the original input plus its raw marker. Since the Huffman
/// candidate is built before this choice is made, it can exceed `max_len` and
/// return an error even when the raw frame would fit.
///
/// A zero `max_len` always returns [`EncodeError::OutputTooLarge`], including
/// for empty input, because even an empty frame requires one byte.
///
/// # Examples
///
/// ```
/// # use huffman_zandronum::encode;
/// let mut output = Vec::new();
/// let len = encode(b"hello", &mut output, 1024)?;
/// assert!(len > 0);
/// # Ok::<(), huffman_zandronum::EncodeError>(())
/// ```
pub fn encode(input: &[u8], output: &mut Vec<u8>, max_len: usize) -> Result<usize, EncodeError> {
    output.clear();
    if max_len == 0 {
        return Err(EncodeError::OutputTooLarge);
    }

    // Create the output already with the padding byte initialized.
    output.push(0);

    // Counts the bits to write a byte
    let mut bit_pos: u8 = 0;

    // Holds the byte we are building
    let mut byte: u8 = 0;

    for &input_byte in input {
        let (code, bit_count) = HUFFMAN_ZANDRONUM.codes[input_byte as usize];

        // Walk the byte from left to right (MSB)
        for cursor in (0..bit_count).rev() {
            // Take the bit, cast u8 as we mask it.
            let bit = ((code >> cursor) & 1) as u8;

            // Active the bit in the position and increment the counter.
            byte |= bit << (7 - bit_pos);
            bit_pos += 1;

            if bit_pos == 8 {
                // We completed the byte
                if output.len() < max_len {
                    output.push(if REVERSE_BITS {
                        byte.reverse_bits()
                    } else {
                        byte
                    });
                } else {
                    return Err(EncodeError::OutputTooLarge);
                }

                bit_pos = 0;
                byte = 0;
            }
        }
    }

    // If the last bit is not complete write the padding and push.
    if bit_pos != 0 {
        output[0] = 8 - bit_pos;
        if output.len() < max_len {
            output.push(if REVERSE_BITS {
                byte.reverse_bits()
            } else {
                byte
            });
        } else {
            return Err(EncodeError::OutputTooLarge);
        }
    }

    if EXPANDABLE {
        return Ok(output.len());
    }

    if output.len() > input.len() + 1 && input.len() + 1 <= max_len {
        output.clear();
        output.push(0xFF);
        output.extend_from_slice(input);
    }

    Ok(output.len())
}

// Golden vectors in this module were produced by the original C++ codec.
#[cfg(test)]
mod tests {
    use super::{DecodeError, EncodeError, decode, encode};

    const MAX_OUTPUT_LEN: usize = 1024;

    struct GoldenVector {
        name: &'static str,
        plaintext: &'static [u8],
        encoded: &'static [u8],
    }

    const GOLDEN_VECTORS: &[GoldenVector] = &[
        GoldenVector {
            name: "empty input",
            plaintext: &[],
            encoded: &[0x00],
        },
        GoldenVector {
            name: "hello uses raw fallback",
            plaintext: b"hello",
            encoded: &[0xff, b'h', b'e', b'l', b'l', b'o'],
        },
        GoldenVector {
            name: "single ASCII byte uses raw fallback",
            plaintext: b"a",
            encoded: &[0xff, b'a'],
        },
        GoldenVector {
            name: "zero byte",
            plaintext: &[0x00],
            encoded: &[0x05, 0x02],
        },
        GoldenVector {
            name: "control byte",
            plaintext: &[0x04],
            encoded: &[0x00, 0xd9],
        },
        GoldenVector {
            name: "space",
            plaintext: b" ",
            encoded: &[0x04, 0x0e],
        },
        GoldenVector {
            name: "pangram",
            plaintext: b"The quick brown fox jumps over the lazy dog",
            encoded: &[
                0x01, 0x2f, 0x66, 0x23, 0x39, 0xb6, 0x31, 0xbf, 0x79, 0x32, 0xf8, 0xe2, 0xc6, 0x51,
                0x53, 0x91, 0xf2, 0xb6, 0xa3, 0x8e, 0x78, 0xc4, 0x31, 0x1b, 0x5a, 0x8d, 0xf5, 0x51,
                0xf8, 0x23, 0x35, 0xfa, 0x5a, 0xb3, 0x91, 0x7c, 0x44, 0xd3, 0x69, 0x9b, 0xbe, 0x71,
                0x94, 0x39,
            ],
        },
        GoldenVector {
            name: "real request",
            plaintext: &[
                0xc7, 0x00, 0x00, 0x00, 0xa0, 0x00, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
            ],
            encoded: &[0x03, 0xb8, 0x49, 0x9a, 0x29, 0x8d, 0x24, 0x09],
        },
        GoldenVector {
            name: "real response",
            plaintext: &[
                0x77, 0x5d, 0x56, 0x00, 0x00, 0x00, 0x00, 0x00, 0x33, 0x2e, 0x33, 0x2d, 0x61, 0x6c,
                0x70, 0x68, 0x61, 0x2d, 0x72, 0x32, 0x36, 0x30, 0x39, 0x32, 0x35, 0x2d, 0x32, 0x33,
                0x35, 0x32, 0x20, 0x6f, 0x6e, 0x20, 0x44, 0x61, 0x72, 0x77, 0x69, 0x6e, 0x20, 0x32,
                0x37, 0x2e, 0x30, 0x2e, 0x30, 0x00, 0xa0, 0x00, 0x08, 0x00, 0x20, 0x00, 0x00, 0x00,
                0x00,
            ],
            encoded: &[
                0x04, 0x53, 0x13, 0xf1, 0x92, 0x24, 0xad, 0xb6, 0xac, 0x36, 0x37, 0x1d, 0xd1, 0xca,
                0x6c, 0xd3, 0xe6, 0xc6, 0xc9, 0xb8, 0x9e, 0x80, 0xc9, 0x9f, 0xe6, 0xc9, 0xd5, 0x9f,
                0x49, 0x1f, 0x85, 0x94, 0x73, 0x6b, 0xda, 0x38, 0xf5, 0x17, 0x29, 0x9f, 0x7c, 0x6e,
                0xe9, 0x69, 0xe9, 0xd1, 0x4c, 0x69, 0xe4, 0x92, 0x04,
            ],
        },
        GoldenVector {
            name: "repeated zero bytes",
            plaintext: &[0x00; 16],
            encoded: &[0x00, 0x92, 0x24, 0x49, 0x92, 0x24, 0x49],
        },
        GoldenVector {
            name: "mixed common bytes",
            plaintext: &[
                0x00, 0x00, 0x00, 0x00, 0x01, 0x02, 0x03, 0x04, 0x20, 0x20, 0x20, 0x20,
            ],
            encoded: &[0x00, 0x92, 0xb4, 0x77, 0x22, 0xd9, 0xee, 0xee],
        },
        GoldenVector {
            name: "mixed binary bytes",
            plaintext: &[0x00, 0x01, 0x02, 0x03, 0x7f, 0xff, 0x20, 0x41],
            encoded: &[0x05, 0xda, 0x3b, 0x91, 0xdc, 0x60, 0xe7, 0x04],
        },
        GoldenVector {
            name: "repeated x uses raw fallback",
            plaintext: &[b'x'; 32],
            encoded: &[
                0xff, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78,
                0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78, 0x78,
                0x78, 0x78, 0x78, 0x78, 0x78,
            ],
        },
        GoldenVector {
            name: "input containing the raw marker",
            plaintext: &[0xff, 0x00, 0x01, 0x02],
            encoded: &[0x05, 0x83, 0xb5, 0x77, 0x02],
        },
    ];

    #[test]
    fn encode_matches_zandronum_golden_vectors() {
        let mut output = Vec::new();

        for case in GOLDEN_VECTORS {
            let len = encode(case.plaintext, &mut output, MAX_OUTPUT_LEN)
                .unwrap_or_else(|error| panic!("{}: {error:?}", case.name));

            assert_eq!(len, case.encoded.len(), "{}: encoded length", case.name);
            assert_eq!(&output[..len], case.encoded, "{}", case.name);
        }
    }

    #[test]
    fn decode_matches_zandronum_golden_vectors() {
        let mut output = Vec::new();

        for case in GOLDEN_VECTORS {
            let len = decode(case.encoded, &mut output, MAX_OUTPUT_LEN)
                .unwrap_or_else(|error| panic!("{}: {error:?}", case.name));

            assert_eq!(len, case.plaintext.len(), "{}: decoded length", case.name);
            assert_eq!(&output[..len], case.plaintext, "{}", case.name);
        }
    }

    #[test]
    fn encode_and_decode_round_trip_golden_inputs() {
        let mut encoded = Vec::new();
        let mut decoded = Vec::new();

        for case in GOLDEN_VECTORS {
            let encoded_len = encode(case.plaintext, &mut encoded, MAX_OUTPUT_LEN)
                .unwrap_or_else(|error| panic!("{}: encode: {error:?}", case.name));
            let decoded_len = decode(&encoded[..encoded_len], &mut decoded, MAX_OUTPUT_LEN)
                .unwrap_or_else(|error| panic!("{}: decode: {error:?}", case.name));

            assert_eq!(&decoded[..decoded_len], case.plaintext, "{}", case.name);
        }
    }

    #[test]
    fn decode_does_not_consume_padding_bits() {
        let mut output = Vec::new();
        let len = decode(&[0x05, 0x0e], &mut output, MAX_OUTPUT_LEN).unwrap();

        assert_eq!(len, 0);
        assert!(output.is_empty());

        let len = decode(&[0x04, 0x0e], &mut output, MAX_OUTPUT_LEN).unwrap();
        assert_eq!(&output[..len], b" ");
    }

    #[test]
    fn output_limits_allow_exact_fit_and_reject_overflow() {
        let mut output = Vec::new();

        let len = encode(&[], &mut output, 1).unwrap();
        assert_eq!(&output[..len], &[0x00]);
        assert!(matches!(
            encode(&[], &mut output, 0),
            Err(EncodeError::OutputTooLarge)
        ));

        let len = decode(&[0xff, b'a'], &mut output, 1).unwrap();
        assert_eq!(&output[..len], b"a");
        assert!(matches!(
            decode(&[0xff, b'a'], &mut output, 0),
            Err(DecodeError::OutputTooLarge)
        ));

        assert!(matches!(
            decode(&[0x05, 0x02], &mut output, 0),
            Err(DecodeError::OutputTooLarge)
        ));
    }
}
