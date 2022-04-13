from bitstring import BitStream

parity_check_mat = BitStream("0b"+
    "1000000000"+
    "0100000000"+
    "0010000000"+
    "0001000000"+
    "0000100000"+
    "0000010000"+
    "0000001000"+
    "0000000100"+
    "0000000010"+
    "0000000001"+
    "1011011100"+
    "0101101110"+
    "0010110111"+
    "1010000111"+
    "1110011111"+
    "1100010011"+
    "1101010101"+
    "1101110110"+
    "0110111011"+
    "1000000001"+
    "1111011100"+
    "0111101110"+
    "0011110111"+
    "1010100111"+
    "1110001111"+
    "1100011011"
)

syndrome_vecs = BitStream("0b"+
    "1111011000"+
    "1111010100"+
    "1001011100"+
    "1111001100"+
    "1001011000"
)

program_type = [
    "No programme type or undefined",
    "News",
    "Information",
    "Sports",
    "Talk",
    "Rock",
    "Classic rock",
    "Adult hits",
    "Soft rock",
    "Top 40",
    "Country",
    "Oldies",
    "Soft music",
    "Nostalgia",
    "Jazz",
    "Classical",
    "Rhythm and blues",
    "Soft rhythm and blues",
    "Language",
    "Religious music",
    "Religious talk",
    "Personality",
    "Public",
    "College",
    "Spanish Talk",
    "Spanish Music",
    "Hip hop",
    "Unassigned",
    "Unassigned",
    "Weather",
    "Emergency test",
    "Emergency"
]

def sync_start(stream):
    idx = 0
    # print(len(stream))
    # print(stream[idx:idx+26])
    offset_type = check_syndrom(stream[idx:idx+26])
    while offset_type != "A" and (idx+26) < len(stream):
        idx += 1
        offset_type = check_syndrom(stream[idx:idx+26])
        # print(offset_type)
    return idx, offset_type

def check_syndrom(msg_block):
    syndrom_vec = BitStream(10)
    for col in range(10):
        parity_vec = BitStream(26)
        for row in range(26):
            parity_vec[row] = '0b1' if parity_check_mat[row*10+col] else '0b0'
        # print(parity_vec)
        # print(len(msg_block))
        prod = msg_block & parity_vec  # element-wise multiplication (AND)
        syndrom_vec[col] = '0b1' if prod.count(1) % 2 else '0b0'  # sum (XOR) of products
    # print(syndrom_vec)

    found = syndrome_vecs.find(syndrom_vec)  # return tuple of index if found, else empty tuple
    if found:
        idx = found[0]
        if idx == 0:
            return "A"
        elif idx == 10:
            return "B"
        elif idx == 20:
            return "Cp"
        elif idx == 30:
            return "C"
        elif idx == 40:
            return "D"
    return -1

def parse_msg(msg_block, block_type):
    if block_type == "A":
        # print(msg_block)
        PI = msg_block[:16]
        # print("PI", len(PI))
        print(f"PI code: {PI.hex}")
    elif block_type == "B":
        group_type = msg_block[:4]
        PTY = program_type[msg_block[6:11].uint]
        DI_seg = msg_block[14:16]
        return group_type, PTY, DI_seg
    elif block_type == "D":
        return chr(msg_block[:8].uint), chr(msg_block[8:16].uint)


