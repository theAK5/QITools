class Qdensitymatrix:

    def __init__(self, rho):
        self.rho = rho
        self.dim = len(rho)
    
    def neg(self,bitstring):
        """
        Function that computes the negativity of a state (density matrix) according to the bitstring. The bitstring is a
        list that determines the required bipartition of the state, and it should consist of only 1 and 0. 

        For example (0,1,0,1,1,0) will compute the negativity when the state is partitioned by grouping qubits 0,2,5 and
        qubits 1,3,4 for a six qubit state.

        Params:
            rho (numpy array): Input density matrix 
            bitstring (list): List that determines the bipartition of the state

        Returns:
            Ng: The negativity of input rho with respect to bipartition specified using bitstring
        """

        M = np.zeros((self.dim,self.dim),dtype = complex)
        for i in range(self.dim):
            for j in range(self.dim):
                a,b = swap(i,j,bitstring)
                M[i][j] = np.copy(rho[a][b])
        
        
        w = np.linalg.eigvalsh(M)
        Ng = np.sum([np.abs(i) for i in w if i<0])
        return(Ng)

    def logneg(self,bitstring):
        """
        Function that computes the logarithmic negativity of the input state rho with respect to bipartition specified b
        bitstring

        Params:
            rho (numpy array): Input density matrix 
            bitstring (list): List that determines the bipartition of the state

        Returns:
            The negativity of input rho with respect to bipartition specified using bitstring
        """
        
        N = neg(self.rho, bitstring) 
        return(np.log2(2*N+1))

    def swap(i1,i2,bitstring):
        
        #Function to help neg function
        bini1 = np.zeros(len(bitstring), dtype = int)
        bini2 = np.zeros(len(bitstring), dtype = int)
        a = i1
        b = i2
        j = k = 0

        while(a>0):
            rem = (a%2)
            a = int(a/2)
            bini1[j] = rem
            j = j + 1
        
        while(b>0):
            rem = (b%2)
            b = int(b/2)
            bini2[k] = rem
            k = k + 1
        
        for i in bitstring[::-1]: 
            if i==1:
                bini1[i],bini2[i] = bini2[i],bini1[i]
        

        deci1 = 0 
        deci2 = 0
        
        for i in range(len(bitstring)):
            deci1 = deci1 + bini1[i]*2**i
            deci2 = deci2 + bini2[i]*2**i
        return(deci1,deci2)

    def parttrace(self,bitstring,sig):
        """
        Function that partially traces a given density matrix according to bitstring and sig. bitstring is a list consisting of 0 and 1. Here 1 denotes the positions of the qubits (qudits) to be traced out. sig is a list of the dimensions of the input state.

        For example: (0,1,0,0,1) will trace out 2nd and 5th qubit/qudit. sig = (2,2,2,2,2) will specify that all qudits are of dimension 2.

        Parameters:
            bitstring(list): List consisting the qubit numbers that are to be traced out. 
            sig(list): List consisting of the dimension of every qubit.

        Returns:
            rho_red (Qdensitymatrix object): The reduced density matrix

        """
        dim = len(self.rho)
        bits = len(bitstring)

        indices = []
        for i in range(bits):
            if(bitstring[i]==1):
                indices.append(i)

        args = []
        for d in sig:
            aux = []
            for i in range(d):
                aux.append(i)
            args.append(aux)
        
        map1 = []
        for combinations in itertools.product(*args):
            map1.append(combinations)

        sig_red = []
        d_reduced = 1
        for i in range(bits):
            if(bitstring[i]!=1):
                d_reduced = d_reduced*sig[i]
                sig_red.append(sig[i])
        args = []
        for d in sig_red:
            aux = []
            for i in range(d):
                aux.append(i)
            args.append(aux)

        map2 = []
        for combinations in itertools.product(*args):
            map2.append(combinations)
        rho_red = np.zeros((d_reduced,d_reduced),dtype = complex)
        for i in range(dim):
            for j in range(dim):
                match = True
                stringi = map1[i] 
                stringj = map1[j]
                strired = list(stringi)
                strjred = list(stringj)
                
                for index in indices:
                    if(stringi[index]!=stringj[index]):
                        match = False

                if(match):
                    for index in sorted(indices,reverse = True):
                        del strired[index]
                        del strjred[index]
                    
                    ir = map2.index(tuple(strired))
                    jr = map2.index(tuple(strjred))
                    rho_red[ir][jr] = rho_red[ir][jr] + self.rho[i][j]
        return(Qdensitymatrix(rho_red))
