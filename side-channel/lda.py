import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import random

def LDA(LDA_traces, LDA_training_coefficient_intermediate_values, coefficient_traces, LDA_parameters, coefficient_intermediate_values=None):
    inv_traces = coefficient_traces[0]
    loop_traces = coefficient_traces[1]
    nb_decapsulation_traces_inv = inv_traces.shape[0]
    nb_decapsulation_traces_loop = loop_traces.shape[0]
    nb_alpha_per_trace = inv_traces.shape[1]
    nb_loop_per_alpha = loop_traces.shape[2]
    inv_traces = inv_traces.reshape(inv_traces.shape[0]*inv_traces.shape[1],inv_traces.shape[-1])
    loop_traces = loop_traces.reshape(loop_traces.shape[0]*loop_traces.shape[1]*loop_traces.shape[2],loop_traces.shape[-1])
    # Testing
    if coefficient_intermediate_values == None:
        lda_inv = np.load(LDA_parameters, allow_pickle=True)
        lda_inv = lda_inv['lda']
        lda_loop = lda_inv[1]
        lda_inv = lda_inv[0]

        inv_traces = lda_inv.transform(inv_traces)
        loop_traces = lda_loop.transform(loop_traces)

        inv_traces = inv_traces.reshape(nb_decapsulation_traces_inv, nb_alpha_per_trace, inv_traces.shape[-1])
        loop_traces = loop_traces.reshape(nb_decapsulation_traces_loop, nb_alpha_per_trace, nb_loop_per_alpha, loop_traces.shape[-1])
        print(inv_traces.shape)
        print(loop_traces.shape)
        
        np.savez_compressed(LDA_test_traces, inv_traces=inv_traces, loop_traces=loop_traces)        
    # Training 
    else:
        inv_intermediate_values = coefficient_intermediate_values[0].reshape(coefficient_intermediate_values[0].shape[0]*coefficient_intermediate_values[0].shape[1])
        loop_intermediate_values = coefficient_intermediate_values[1].reshape(coefficient_intermediate_values[1].shape[0]*coefficient_intermediate_values[1].shape[1]*coefficient_intermediate_values[1].shape[2])
        print(inv_intermediate_values.shape)
        print(loop_intermediate_values.shape)
        inv_hw_intermediate_values = np.array([bin(intermediate_value).count('1') for intermediate_value in inv_intermediate_values])
        print(inv_hw_intermediate_values.shape)
        loop_hw_intermediate_values = np.array([bin(intermediate_value).count('1') for intermediate_value in loop_intermediate_values])
        print(loop_hw_intermediate_values.shape)

        lda_inv = LinearDiscriminantAnalysis()
        lda_inv.fit(inv_traces, inv_hw_intermediate_values)
        inv_traces = lda_inv.transform(inv_traces)

        lda_loop = LinearDiscriminantAnalysis()
        lda_loop.fit(loop_traces, loop_hw_intermediate_values)
        loop_traces = lda_loop.transform(loop_traces)
        
        print(inv_traces.shape)
        print(loop_traces.shape)
        
        inv_traces = inv_traces.reshape(nb_decapsulation_traces_inv, nb_alpha_per_trace, inv_traces.shape[-1])
        loop_traces = loop_traces.reshape(nb_decapsulation_traces_loop, nb_alpha_per_trace, nb_loop_per_alpha, loop_traces.shape[-1])
        print(inv_traces.shape)
        print(loop_traces.shape)
        
        trained_lda = np.zeros((2), dtype=type(LinearDiscriminantAnalysis()))
        trained_lda[0] = lda_inv
        trained_lda[1] = lda_loop
        np.savez_compressed(LDA_parameters, lda=trained_lda)
        
        np.savez_compressed(LDA_training_traces, inv_traces=inv_traces, loop_traces=loop_traces)
        np.savez_compressed(LDA_training_coefficient_intermediate_values, inv_intermediate_values=coefficient_intermediate_values[0], loop_intermediate_values=coefficient_intermediate_values[1])


if __name__ == '__main__':
    import argparse
 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-t','--training', action='store_true', dest='training', default=False)
    
    parser.add_argument('--me','--mceliece', action='store', dest='mceliece', required=True,
                        choices=['51220', '348864', '8192128'],help='McEliece security paramaters')

    args = parser.parse_args()
    
    coefficient_traces = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_traces_'+str(args.mceliece)+'.npz'
    coefficient_intermediate_values = 'results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'.npz'

    LDA_parameters = 'results/matrices/'+str(args.mceliece)+'/attack/LDA_parameters_'+str(args.mceliece)+'.npz'
    
    LDA_training_traces = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_traces_'+str(args.mceliece)+'.npz'
    LDA_training_coefficient_intermediate_values = 'results/matrices/'+str(args.mceliece)+'/attack/lda_training_coefficient_intermediate_values_'+str(args.mceliece)+'.npz'

    LDA_test_traces = 'results/matrices/'+str(args.mceliece)+'/attack/lda_test_traces_'+str(args.mceliece)+'.npz'

    if args.mceliece == '8192128':
        list_coefficient_traces = list('results/matrices/'+str(args.mceliece)+'/attack/coefficient_traces_'+str(args.mceliece)+'_'+str(nb_file)+'.npz' for nb_file in range(5))
        coefficient_traces = np.load(list_coefficient_traces[0])
        inv_traces = coefficient_traces['inv_traces_matrix']
        loop_traces = coefficient_traces['loop_traces_matrix'][:5]
        loop_traces = np.concatenate((loop_traces, coefficient_traces['loop_traces_matrix'][7:8]), axis=0)
        loop_traces = np.concatenate((loop_traces, coefficient_traces['loop_traces_matrix'][31:32]), axis=0)
        loop_test_traces = coefficient_traces['loop_traces_matrix'][49:51]
        for path_to_matrix in list_coefficient_traces[1:]:
            with np.load(path_to_matrix) as coefficient_traces:
                inv_traces = np.concatenate((inv_traces, coefficient_traces['inv_traces_matrix']), axis=0)
        print(inv_traces.shape)
        print(loop_traces.shape)
        list_coefficient_intermediate_values = list('results/matrices/'+str(args.mceliece)+'/attack/coefficient_intermediates_values_'+str(args.mceliece)+'_'+str(nb_file)+'.npz' for nb_file in range(5))
        inv_intermediate_values = np.load( list_coefficient_intermediate_values[0])
        loop_intermediate_values = inv_intermediate_values['loop_intermediate_values'][:5]
        loop_intermediate_values = np.concatenate((loop_intermediate_values, inv_intermediate_values['loop_intermediate_values'][7:8]), axis=0)
        loop_intermediate_values = np.concatenate((loop_intermediate_values, inv_intermediate_values['loop_intermediate_values'][31:32]), axis=0)
        inv_intermediate_values = inv_intermediate_values['inv_intermediate_values']
        for path_to_matrix in list_coefficient_intermediate_values[1:]:
            with np.load(path_to_matrix) as coefficient_intermediate_values:
                inv_intermediate_values = np.concatenate((inv_intermediate_values, coefficient_intermediate_values['inv_intermediate_values']), axis=0)
        print(inv_intermediate_values.shape)
        print(loop_intermediate_values.shape)
    else:
        coefficient_traces = np.load(coefficient_traces)
        inv_traces = coefficient_traces['inv_traces_matrix']
        loop_traces = coefficient_traces['loop_traces_matrix']
        inv_intermediate_values = np.load(coefficient_intermediate_values)
        loop_intermediate_values = inv_intermediate_values['loop_intermediate_values']
        inv_intermediate_values = inv_intermediate_values['inv_intermediate_values']

    if args.mceliece == '51220':
        nb_training_traces_inv = 200
        nb_training_traces_loop = 200
        nb_test_traces_inv = nb_training_traces
        nb_test_traces_loop = nb_testing_traces_loop
    if args.mceliece == '348864':
        nb_training_traces_inv = 500
        nb_training_traces_loop = 300
        nb_test_traces_inv = 500
        nb_test_traces_loop = 500
    if args.mceliece == '8192128':
        nb_training_traces_inv = 700
        nb_test_traces_inv = nb_training_traces_inv
    
    if args.mceliece == '8192128':
        if args.training: LDA(LDA_training_traces, LDA_training_coefficient_intermediate_values, (inv_traces[:nb_training_traces_inv], loop_traces), LDA_parameters, (inv_intermediate_values[:nb_training_traces_inv], loop_intermediate_values))
        else: LDA(LDA_test_traces, LDA_training_coefficient_intermediate_values, (inv_traces[nb_test_traces_inv:], loop_test_traces), LDA_parameters)
    else :
        if args.training: LDA(LDA_training_traces, LDA_training_coefficient_intermediate_values, (inv_traces[:nb_training_traces_inv], loop_traces[:nb_training_traces_loop]), LDA_parameters, (inv_intermediate_values[:nb_training_traces_inv], loop_intermediate_values[:nb_training_traces_loop]))
        else: LDA(LDA_test_traces, LDA_training_coefficient_intermediate_values, (inv_traces[nb_test_traces_inv:], loop_traces[nb_test_traces_loop:]), LDA_parameters)