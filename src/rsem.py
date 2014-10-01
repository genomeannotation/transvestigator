def read_rsem(io_buffer):
    data = []
    columns = io_buffer.readline()
    for line in io_buffer:
        cols = line.strip().split("\t")
        data.append((cols[0], cols[1], int(cols[2]), float(cols[3]), float(cols[4]), float(cols[5]), float(cols[6]), float(cols[7])))
    return data
