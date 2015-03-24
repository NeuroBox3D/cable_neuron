NEURON {
    POINT_PROCESS rel_exp
    RANGE tau1, tau2
    
    RANGE T
    GLOBAL total
}

UNITS {    
    (mM) = (milli/liter)
}

PARAMETER {
    tau1=.1 (ms) <1e-9,1e9>
    tau2 = 10 (ms) <1e-9,1e9>
}

ASSIGNED {
    T (mM) 
    factor
    total (mM)
}

STATE {
    A (mM)
    B (mM)
}

INITIAL {
    LOCAL tp
    total = 0
    if (tau1/tau2 > .9999) {
	tau1 = .9999*tau2
    }
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    T = B - A
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
    state_discontinuity(A, A + weight*factor)
    state_discontinuity(B, B + weight*factor)
    total = total+weight
}