import ast

def read_config(io_buffer):
    config = {}
    for line in io_buffer:
        if line.startswith('#') or line.strip() == "":
            continue
        key, val = line.split("=")
        key = key.strip()
        val = val.strip()
        config[key] = ast.literal_eval(val)
    return config
