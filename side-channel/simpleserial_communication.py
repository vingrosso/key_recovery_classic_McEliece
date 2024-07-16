from math import *
import sys
from numpy import uint8

def set_parameter_on_target(target, simpleserial_cmd, parameter_value,nb_values_max=32):
    #simpleserial uses bytearrays
    if type(parameter_value[0]) == str: parameter_value = bytes.fromhex(parameter_value)
    elif type(parameter_value[0]) == uint8: parameter_value = bytes(parameter_value)
    offset = 0
    nb_values = 0
    msg = bytearray([offset,nb_values]+[0]*nb_values_max)
    while (offset * nb_values_max) + nb_values < len(parameter_value):
        msg[2 + nb_values]=parameter_value[(offset * nb_values_max) + nb_values]
        nb_values += 1
        msg[1] = nb_values
        if nb_values == nb_values_max:
            target.simpleserial_write(simpleserial_cmd, msg)
            msg_received = target.simpleserial_read(simpleserial_cmd, nb_values_max, timeout=0)
            if msg[2:] != msg_received:
                print(msg)
                print(msg_received)
                return False
            offset += 1
            nb_values = 0
            msg = bytearray([offset,nb_values]+[0]*nb_values_max)
    if nb_values != 0 :
        target.simpleserial_write(simpleserial_cmd, msg)
        msg_received = target.simpleserial_read(simpleserial_cmd, nb_values_max, timeout=0)
        if msg[2:] != msg_received:
            print(msg)
            print(msg_received)
            return False
    return True

def set_parameter_on_target_offset_on_two_bytes(target, simpleserial_cmd, parameter_value,nb_values_max=30):
    #simpleserial uses bytearrays
    if type(parameter_value[0]) == str: parameter_value = bytes.fromhex(parameter_value)
    elif type(parameter_value[0]) == uint8: parameter_value = bytes(parameter_value)
    offset = 0
    nb_values = 0
    msg = bytearray([offset.to_bytes(2, 'little')[0], offset.to_bytes(2, 'little')[1], nb_values]+[0]*nb_values_max)
    while (offset * nb_values_max) + nb_values < len(parameter_value):
        msg[3 + nb_values]=parameter_value[(offset * nb_values_max) + nb_values]
        nb_values += 1
        msg[2] = nb_values
        if nb_values == nb_values_max:
            target.simpleserial_write(simpleserial_cmd, msg)
            msg_received = target.simpleserial_read(simpleserial_cmd, nb_values_max, timeout=0)
            if msg[3:] != msg_received:
                print(msg)
                print(msg_received)
                return False
            offset += 1
            nb_values = 0
            msg = bytearray([offset.to_bytes(2, 'little')[0], offset.to_bytes(2, 'little')[1], nb_values]+[0]*nb_values_max)
    if nb_values != 0 :
        target.simpleserial_write(simpleserial_cmd, msg)
        msg_received = target.simpleserial_read(simpleserial_cmd, nb_values_max, timeout=0)
        if msg[3:] != msg_received:
            print(msg)
            print(msg_received)
            return False
    return True

def support_generation(target, size_of_L_in_bytes, packet_size_in_bytes=32):
    print("Generation of the support L")
    nb_packets = ceil(size_of_L_in_bytes/packet_size_in_bytes)
    target.simpleserial_write('l', bytearray([nb_packets]))
    for i in range(nb_packets):
        print(target.simpleserial_read('l', packet_size_in_bytes, timeout=0, ack=False))
    return True
