import numpy as np

class qcircuit:
    def __init__(self, qubits):
        #initialising the number of qubits and the gates as identity
        self.qubits = qubits
        self.unitary = np.eye(2**qubits)
        self.index = 2**qubits 
    
    def H(self,n):
        #Hadamard gate on line n
        Had = (1/np.sqrt(2))*np.array([[1,1],[1,-1]])
        single_gates = [Had if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
    
    def T(self,n):
        #Hadamard circuit on line n
        T = np.array([[1,0],[0,np.cos(np.pi/4)+1j*np.sin(np.pi/4)]])
        single_gates = [T if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
    
    def X(self,n):
        X = np.array([[0,1],[1,0]])
        single_gates = [X if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
   
    def Y(self,n):
        #Y gate on line n
        Y = np.array([[0,-1j],[1j,0]])
        single_gates = [Y if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
    
    def Z(self,n):
        #X gate on line n
        Z = np.array([[1,0],[0,-1]])
        single_gates = [Z if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
    
    def CNOT(self,n1,n2):
        #CNOT gate on line n1 and n2
        X = np.array([[0,1],[1,0]])
        P0 = np.array([[0,1],[0,0]])
        P1 = np.array([[0,0],[1,0]])
        CN0 = 1
        CN1 = 1 
        
        for i in range(self.qubits):
            if(i==n1):
                CN0 = np.kron(CN0,P0)
                CN1 = np.kron(CN1,P1)
            elif(i==n2):
                CN0 = np.kron(CN0,X)
                CN1 = np.kron(CN1,np.eye(2))
            else:
                CN0 = np.kron(CN0,np.eye(2))
                CN1 = np.kron(CN1,np.eye(2))
        
        Gate = CN0 + CN1
        self.unitary = np.matmul(self.unitary,Gate)
    
    def S(self,n):
        #S gate line n
        S = np.array([[1,0],[0,1j]])
        single_gates = [S if i==n else np.eye(2) for i in range(self.qubits)]
        Gate = 1
        for g in single_gates:
            Gate = np.kron(Gate,g)
        self.unitary = np.matmul(self.unitary,Gate)
    
    def Gates(self):
        #returns the gate set
        return(self.unitary)
