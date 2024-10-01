import sys
import os
import random

def generate_lookup_table(filename='lookup_table.txt'):
    random.seed(253)
    #fixing the seed to generate a repeatable lookup table of seeds
    
    seed_range = (1, 1000000000) 
    #2**32 is roughly 10**9.68
    lookup_table = {i: random.randint(*seed_range) for i in range(1, 50001)}
    with open(filename, 'w') as file:
        for serial_num, seed in lookup_table.items():
            file.write(f"{serial_num},{seed}\n")
    #os.chmod(filename, 0o444)
    return lookup_table

def read_seed(serial_number, filename='lookup_table.txt'):
    with open(filename, 'r') as file:
        for line in file:
            serial, seed = line.strip().split(',')
            if int(serial) == serial_number:
                return int(seed)
    return None

def get_seed(serial_number):
    lookup_table = generate_lookup_table()
    seed = lookup_table[serial_number]
    return seed
    
if __name__ == "__main__":
    if len(sys.argv) == 2:
        serial_number = int(sys.argv[1])
        seed = get_seed(serial_number)
        if seed is not None:
            print(serial_number, seed)
        else:
            print("Serial number not found")
    else:
        print("Usage: python lookup_script.py <serial_number>")

