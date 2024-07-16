import numpy as np
import chipwhisperer as cw
import subprocess
from math import *
import sys
import argparse
from simpleserial_communication import * 
from tqdm import tqdm

def traces_acquisition(traces_matrix_output, parameters_input, chipwhisperer, stream_mode, compilation, leakage_assessment, mceliece, check_synd):
    scope = cw.scope()
    target = cw.target(scope)
    scope.default_setup()

    scope.clock.adc_src="clkgen_x1"
    scope.gain.db = 32.3
    
    scope.adc.decimate = 1

    if chipwhisperer == 'pro':
        if stream_mode:
            scope.adc.stream_mode = True
            if leakage_assessment:
                scope.adc.samples = 60_000
            else:
                if mceliece == "51220":
                    scope.adc.samples = 5_250_000
                    scope.adc.offset = 0
                elif mceliece == "348864":
                    scope.adc.samples = 12_500_000
                    scope.adc.offset = 25800
                elif mceliece == "8192128":
                    scope.adc.samples = 12_500_000
                    scope.adc.offset = 30000
        else:
            scope.adc.samples = 96_000
    elif chipwhisperer == 'lite':
        scope.adc.samples = 24_400

    if compilation:
        if stream_mode:
            if check_synd:
                subprocess.call(["make","allquick", "PLATFORM=CW308_STM32F3", "CRYPTO_TARGET=NONE", "OPT=s", "CFLAGS_LAST+='-DCRYPTO_NAMESPACE(x)=x' '-D_CRYPTO_NAMESPACE(x)=_##x' '-DSTREAM_MODE' '-DCHECK_SYND'"],
                                cwd="embedded_code/"+mceliece)
            else:
                subprocess.call(["make","allquick", "PLATFORM=CW308_STM32F3", "CRYPTO_TARGET=NONE", "OPT=s", "CFLAGS_LAST+='-DCRYPTO_NAMESPACE(x)=x' '-D_CRYPTO_NAMESPACE(x)=_##x' '-DSTREAM_MODE'"],
                                cwd="embedded_code/"+mceliece)
        else:
            subprocess.call(["make","allquick", "PLATFORM=CW308_STM32F3", "CRYPTO_TARGET=NONE", "OPT=s", "CFLAGS_LAST+='-DCRYPTO_NAMESPACE(x)=x' '-D_CRYPTO_NAMESPACE(x)=_##x'"],
                            cwd="embedded_code/"+mceliece)

    #Uploading the program to the target

    cw.program_target(scope, cw.programmers.STM32FProgrammer, "embedded_code/"+mceliece+"/syndrome_computation_"+mceliece+"-CW308_STM32F3.hex")

    #Cleaning project

    if compilation:
        subprocess.call(["make", "PLATFORM=CW308_STM32F3", "clean"],
                        cwd="embedded_code/"+mceliece)

    #Preparation of the capture

    with np.load(parameters_input) as parameters_matrices:
        if leakage_assessment : 
            nb_traces = len(parameters_matrices["ct_matrix"])
            nb_max_traces_per_matrix = 1000
        else : 
            if mceliece == "51220":
                nb_traces = 300
                nb_max_traces_per_matrix = 200
            elif mceliece == "348864":
                nb_traces = 605
                nb_max_traces_per_matrix = 100
            elif mceliece == "8192128":
                nb_traces = 1000
                nb_max_traces_per_matrix = 100
        nb_matrices = 0
        if nb_traces < nb_max_traces_per_matrix : traces_matrix = np.zeros((nb_traces, scope.adc.samples), dtype="float")
        else: traces_matrix = np.zeros((nb_max_traces_per_matrix, scope.adc.samples), dtype="float")
        print(nb_traces)
        print(traces_matrix.shape)
        nb_row = 0
        ct_set = False
        g_set = False
        L_set = False
        sys_t = len(parameters_matrices["g_matrix"][0])//2-1
        for index_trace in tqdm(range(nb_traces)):
            if stream_mode:
                if not np.any(parameters_matrices["ct_matrix"][index_trace]) : sys.exit("Error : There is no ciphered text ct to send")
                if mceliece == "51220" :
                    if not set_parameter_on_target(target,  'r', parameters_matrices["ct_matrix"][index_trace]) : sys.exit("Error : The ciphered text ct has not been sent to the target")
                elif mceliece == "348864" :
                    if not set_parameter_on_target(target,  'r', parameters_matrices["ct_matrix"][index_trace,:13]) : sys.exit("Error : The ciphered text ct has not been sent to the target")
                elif mceliece == "8192128" :
                    if not set_parameter_on_target(target,  'r', parameters_matrices["ct_matrix"][index_trace,:10]) : sys.exit("Error : The ciphered text ct has not been sent to the target")

                if not np.any(parameters_matrices["g_matrix"][index_trace]): sys.exit("Error : There is no goppa polynomial g to send")
                if not set_parameter_on_target(target, 'g', parameters_matrices["g_matrix"][index_trace]) : sys.exit("Error : The goppa polynomial g has not been sent to the target")

                if not np.any(parameters_matrices["L_matrix"][index_trace]) : sys.exit("Error : There is no support L to send")
                if mceliece == "51220" :
                    if not set_parameter_on_target(target, 'L', parameters_matrices["L_matrix"][index_trace]) : sys.exit("Error : The support L has not been sent to the target")
                elif mceliece == "348864" :
                    if not set_parameter_on_target(target, 'L', parameters_matrices["L_matrix"][index_trace,:204]) : sys.exit("Error : The support L has not been sent to the target")                    
                elif mceliece == "8192128":
                    if not set_parameter_on_target_offset_on_two_bytes(target, 'L', parameters_matrices["L_matrix"][index_trace,:154]) : sys.exit("Error : The support L has not been sent to the target")
            else:
                if not np.any(parameters_matrices["ct_matrix"][index_trace]) : sys.exit("Error : There is no ciphered text ct to send")
                if not set_parameter_on_target(target, 'r', [parameters_matrices["ct_matrix"][index_trace][0]],1) : sys.exit("Error : The ciphered text ct has not been sent to the target")

                if not np.any(parameters_matrices["g_matrix"][index_trace]): sys.exit("Error : There is no goppa polynomial g to send")
                if not set_parameter_on_target(target,'g',parameters_matrices["g_matrix"][index_trace]) : sys.exit("Error : The goppa polynomial g has not been sent to the target")

                if not np.any(parameters_matrices["L_matrix"][index_trace]) : sys.exit("Error : There is no support L to send")
                if not set_parameter_on_target(target,'L',[parameters_matrices["L_matrix"][index_trace][0],parameters_matrices["L_matrix"][index_trace][1]],2) : sys.exit("Error : The support L has not been sent to the target")

            scope.arm()
            target.simpleserial_write('s', bytearray())
            scope.capture()

            if check_synd:
                sd_target = bytearray()
                for i in range((sys_t*4)//32):
                    msg_received = target.simpleserial_read('s', 32, timeout=0, ack=False)
                    sd_target += msg_received
                if (sys_t*4)%32:
                    msg_received = target.simpleserial_read('s', 32, timeout=0, ack=False)
                    sd_target += msg_received[:(sys_t*4)%32]
                if not np.any(parameters_matrices["sd_matrix"][index_trace]):
                    sys.exit("Error : No syndrome sd to compare the computed syndrome with")
                sd = bytearray(parameters_matrices["sd_matrix"][index_trace])
                if sd != sd_target :
                    print(sd)
                    print(sd_target)
                    sys.exit("Error : Syndroms do not match")
            else:
                for i in range(1):
                    msg_received = target.simpleserial_read('s', 1, timeout=0, ack=False)
            traces_matrix[nb_row][:] = scope.get_last_trace()
            nb_row += 1                
            if (nb_row >= nb_max_traces_per_matrix) or (index_trace+1 == nb_traces) :
                np.savez_compressed(traces_matrix_output+str(nb_matrices)+".npz", traces_matrix=traces_matrix)
                nb_matrices += 1
                nb_row = 0
                if (index_trace+1 < nb_traces) and (index_trace + 1 + nb_max_traces_per_matrix > nb_traces):
                    traces_matrix = np.zeros((nb_traces-(index_trace+1), scope.adc.samples), dtype="float")
            ct_set = False
            g_set = False
            L_set = False
    scope.dis()
    target.dis()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--cw', action='store', dest='chipwhisperer', required=True,
                        choices=['pro', 'lite'],help='Chipwhisperer used')

    stream_mode = parser.add_mutually_exclusive_group(required=True)
    stream_mode.add_argument('-s','--streaming', dest='stream_mode', action='store_true')
    stream_mode.add_argument('--ns','--no-streaming', dest='stream_mode', action='store_false')

    compilation = parser.add_mutually_exclusive_group(required=True)
    compilation.add_argument('-c','--compilation', dest='compilation', action='store_true')
    compilation.add_argument('--nc','--no-compilation', dest='compilation', action='store_false')     

    leakage_assessment_or_attack = parser.add_mutually_exclusive_group(required=True)
    leakage_assessment_or_attack.add_argument('--la','--leakage_assessment', dest='leakage_assessment', action='store_true')
    leakage_assessment_or_attack.add_argument('-a','--attack', dest='leakage_assessment', action='store_false')

    check_synd = parser.add_mutually_exclusive_group(required=True)
    check_synd.add_argument('--csd','--check_synd', dest='check_synd', action='store_true')
    check_synd.add_argument('--ncsd','--no-check_synd', dest='check_synd', action='store_false')
    
    parser.add_argument('--me','--mceliece', action='store', dest='mceliece', required=True,
                        choices=['51220', '348864', '8192128'],help='McEliece security paramaters')

    args = parser.parse_args()
    
    if args.leakage_assessment:
        parameters_input = "results/matrices/51220/parameters_synd_51220.npz"
        traces_matrix_output = "results/matrices/51220/decapsulation_traces_matrix_51220_leakage_assessment_"
    else:
        parameters_input = "results/matrices/"+args.mceliece+"/parameters_synd_"+args.mceliece+".npz"
        traces_matrix_output = "results/matrices/"+args.mceliece+"/decapsulation_traces_matrix_"+args.mceliece+"_attack_"

    print(traces_matrix_output)
    print(parameters_input)
    traces_acquisition(traces_matrix_output, parameters_input, args.chipwhisperer, args.stream_mode, args.compilation, args.leakage_assessment, args.mceliece, args.check_synd)