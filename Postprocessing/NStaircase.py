###
# Imports
###
import numpy as np
import scipy.optimize as optimize
import scipy.special as sp
import scipy.integrate as integrate
import copy


###
# Properties
###
class PopProp:

    def __init__(self):
        self.InitCellsNum = 1e2
        self.InitCellsPos = 0
        self.InitCellsGen = 0
        self.BirthRate = 1
        self.StressBirthRate = 0
        self.ResistCost = 0
        self.DeathRate = 1e-1
        self.StressDeathRate = 0
        self.MuteUp = 1e-7
        self.MuteDown = 1e-4
        self.StressMuteUp = 0
        self.StressMuteDown = 0
        self.Move = 1e-3
        self.Chemotax = 0.5
        self.StressDepChemotax = False
        self.Chemokin = 0
        self.SwitchUp = 0
        self.SwitchDown = 0
        self.ConsiderSwarm = False
        self.SwarmingMove = 1e-1
        self.SwarmingDensity = 5e4
        self.ConsiderHGT = False
        self.HGTRate = 0


class LatProp:

    def __init__(self):
        self.GridBound = 8
        self.GenBound = 2
        self.PopulationNumber = 1
        self.CarryingCapacity = 1E5
        self.LD = 1


# Main class
class Staircase:

    def __init__(self, l_prop, p_prop):
        #   Save constants
        self.PProp = copy.deepcopy(p_prop)
        self.LProp = copy.deepcopy(l_prop)
        self.Initialize()

    # Function to initialize, or re-initialize the staircase object
    def Initialize(self):
        #   Compute stable wild type
        self.LProp.GenBound = 1
        N_init0 = []
        N_init1 = []
        for pop in range(self.LProp.PopulationNumber):
            N0 = self.LProp.CarryingCapacity * (1 - self.PProp[pop].DeathRate / self.PProp[pop].BirthRate)
            # various initial guess might work for various situations
            N_init0 = N_init0 + [N0 / self.LProp.PopulationNumber for pos in range(self.LProp.GridBound)]
            N_init1 = N_init1 + [N0/self.LProp.PopulationNumber for pos in range(self.LProp.LD)] + [0 for pos in range(self.LProp.GridBound-self.LProp.LD)]
        N_init0 = np.asarray(N_init0)
        N_init1 = np.asarray(N_init1)
        self.N = optimize.fsolve(self.Force, N_init0)
        if self.N[self.LProp.GridBound-1] <= 0:
            self.N = optimize.fsolve(self.Force, N_init1)
        self.LProp.GenBound = 2

        #   Compute mutant birth rate r_m, and mutant net birth rate r_net
        self.r_m = np.zeros((self.LProp.PopulationNumber, self.LProp.GridBound))
        self.r_net = np.zeros((self.LProp.PopulationNumber, self.LProp.GridBound))
        for pop in range(self.LProp.PopulationNumber):
            for pos in range(self.LProp.GridBound):
                    self.r_m[pop][pos] = self.BirthPC(pop, 1, pos, self.N)
                    self.r_net[pop][pos] = self.BirthPC(pop, 1, pos, self.N)-self.DeathPC(pop, 1, pos)

        #   Compute adaptation rate parameters
        self.nDeff = np.zeros(self.LProp.PopulationNumber)
        self.nOeff = np.zeros(self.LProp.PopulationNumber)
        self.nNeff = np.zeros(self.LProp.PopulationNumber)
        self.aD = np.zeros(self.LProp.PopulationNumber)
        self.bD = np.zeros(self.LProp.PopulationNumber)
        self.cD = np.zeros(self.LProp.PopulationNumber)
        self.lambdaO = np.zeros(self.LProp.PopulationNumber)
        self.aN = np.zeros(self.LProp.PopulationNumber)
        self.bN = np.zeros(self.LProp.PopulationNumber)
        self.cN = np.zeros(self.LProp.PopulationNumber)
        D = self.LProp.LD-1
        O = self.LProp.LD
        N = self.LProp.LD+1
        for pop in range(self.LProp.PopulationNumber):
            if D < self.LProp.GridBound:
                self.nDeff[pop] = self.N[pop*self.LProp.GridBound+D]
                self.aD[pop] = self.MuteUpPC(pop,0,D)*self.nDeff[pop]/(self.r_m[pop][D]+self.PProp[pop].HGTRate*self.nDeff[pop]-self.MuteUpPC(pop,0,D))
                self.bD[pop] = (self.MoveRightPC(pop,1,D,self.N)+self.DeathPC(pop,1,D)-self.PProp[pop].HGTRate*self.nDeff[pop]-self.r_m[pop][D]+self.MuteUpPC(pop,0,D)+self.MuteDownPC(pop,1,D))/2
                self.cD[pop] = np.sqrt((self.MoveRightPC(pop,1,D,self.N)+self.DeathPC(pop,1,D)+self.PProp[pop].HGTRate*self.nDeff[pop]+self.r_m[pop][D]-self.MuteUpPC(pop,0,D)+self.MuteDownPC(pop,1,D))**2+4*(self.MuteUpPC(pop,0,D)-self.PProp[pop].HGTRate*self.nDeff[pop]-self.r_m[pop][D])*(self.DeathPC(pop,1,D)+self.MuteDownPC(pop,1,D)))/2
            else:
                self.nDeff[pop] = 0
                self.aD[pop] = 0
                self.bD[pop] = 1
                self.cD[pop] = 1
            if O < self.LProp.GridBound:
                self.nOeff[pop] = self.N[pop*self.LProp.GridBound+O]
                self.lambdaO[pop] = self.MuteUpPC(pop, 0, O) * self.nOeff[pop]
            else:
                self.nOeff[pop] = 0
                self.lambdaO[pop] = 0
            if N < self.LProp.GridBound:
                self.nNeff[pop] = self.N[pop*self.LProp.GridBound+N]
                self.aN[pop] = self.MuteUpPC(pop,0,N)*self.nNeff[pop]/(self.r_m[pop][N]+self.PProp[pop].HGTRate*self.nNeff[pop]-self.MuteUpPC(pop,0,N))
                self.bN[pop] = (self.MoveLeftPC(pop,1,N,self.N)+self.DeathPC(pop,1,N)-self.PProp[pop].HGTRate*self.nNeff[pop]-self.r_m[pop][N]+self.MuteUpPC(pop,0,N)+self.MuteDownPC(pop,1,N))/2
                self.cN[pop] = np.sqrt((self.MoveLeftPC(pop,1,N,self.N)+self.DeathPC(pop,1,N)+self.PProp[pop].HGTRate*self.nNeff[pop]+self.r_m[pop][N]-self.MuteUpPC(pop,0,N)+self.MuteDownPC(pop,1,N))**2+4*(self.MuteUpPC(pop,0,N)-self.PProp[pop].HGTRate*self.nNeff[pop]-self.r_m[pop][N])*(self.DeathPC(pop,1,N)+self.MuteDownPC(pop,1,N)))/2
            else:
                self.nNeff[pop] = 0
                self.aN[pop] = 0
                self.bN[pop] = 1
                self.cN[pop] = 1

    # Compute N space tot
    def NSpaceTot(self, pos, N):
        N_tot = 0
        gen_num = int(len(N)/(self.LProp.PopulationNumber * self.LProp.GridBound))
        for gen in range(gen_num):
            for pop in range(self.LProp.PopulationNumber):
                N_tot += N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return N_tot

    # Compute N pop space tot
    def NPopSpaceTot(self, pos, pop, N):
        N_tot = 0
        gen_num = int(len(N)/(self.LProp.PopulationNumber * self.LProp.GridBound))
        for gen in range(gen_num):
            N_tot += N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return N_tot

    # Compute per capita birth rate
    def BirthPC(self, pop, gen, pos, N):
        rate = 0
        N_pos = self.NSpaceTot(pos, N)
        if N_pos < self.LProp.CarryingCapacity:
            rate = (1 - N_pos / self.LProp.CarryingCapacity)
            if pos < self.LProp.LD + gen:
                rate *= (1 - self.PProp[pop].ResistCost) ** (self.LProp.LD + gen - pos - 1)
                rate *= self.PProp[pop].BirthRate
            else:
                rate *= self.PProp[pop].StressBirthRate
        return rate

    # Compute death rate per capita
    def DeathPC(self, pop, gen, pos):
        rate = self.PProp[pop].DeathRate
        if pos >= self.LProp.LD + gen:
            rate += self.PProp[pop].StressDeathRate
        return rate

    # Compute mutation up per capita
    def MuteUpPC(self, pop, gen, pos):
        rate = self.PProp[pop].MuteUp
        if pos >= self.LProp.LD + gen:
            rate += self.PProp[pop].StressMuteUp
        return rate

    # Compute HGT rate
    def MuteHGTPC(self, pop, gen, gen_donor, pos, N):
        rate = 0
        if self.PProp[pop].ConsiderHGT:
            rate = self.PProp[pop].HGTRate
            rate *= N[gen_donor * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute mutation down per capita
    def MuteDownPC(self, pop, gen, pos):
        rate = self.PProp[pop].MuteDown
        if pos >= self.LProp.LD + gen:
            rate += self.PProp[pop].StressMuteDown
        return rate

    # Compute move right per capita
    def MoveRightPC(self, pop, gen, pos, N):
        rate = self.PProp[pop].Move
        if self.PProp[pop].ConsiderSwarm and self.PProp[pop].SwarmingDensity <= self.NPopSpaceTot(pos, pop, N):    # if swarming at large enough density
            rate = self.PProp[pop].SwarmingMove
        #if pos > self.LProp.LD + gen - 2: # if to be stressed
        #    rate += self.PProp[pop].Chemokin
        if pos > self.LProp.LD + gen - 1:   # if stressed
            rate += self.PProp[pop].Chemokin
            rate *= 2*self.PProp[pop].Chemotax
        elif not self.PProp[pop].StressDepChemotax:    # if not stressed and stress indep. chemotaxis
            rate *= 2*self.PProp[pop].Chemotax
        return rate

    # Compute move left per capita
    def MoveLeftPC(self, pop, gen, pos, N):
        rate = self.PProp[pop].Move
        if self.PProp[pop].ConsiderSwarm and self.PProp[pop].SwarmingDensity <= self.NPopSpaceTot(pos, pop, N):    # if swarming at large enough density
            rate = self.PProp[pop].SwarmingMove
        #if pos > self.LProp.LD + gen: # if to be stressed
        #    rate += self.PProp[pop].Chemokin
        if pos > self.LProp.LD + gen - 1:   # if stressed
            rate += self.PProp[pop].Chemokin
            rate *= 2*(1 - self.PProp[pop].Chemotax)
        elif not self.PProp[pop].StressDepChemotax:    # if not stressed and stress indep. chemotaxis
            rate *= 2*(1 - self.PProp[pop].Chemotax)
        return rate

    # Compute switch up per capita
    def SwitchUpPC(self, pop, gen, pos):
        rate = self.PProp[pop].SwitchUp
        return rate

    # Compute switch down per capita
    def SwitchDownPC(self, pop, gen, pos):
        rate = self.PProp[pop].SwitchDown
        return rate

    # Compute birth rate
    def Birth(self, pop, gen, pos, N):
        rate = self.BirthPC(pop, gen, pos, N)
        rate *= N[gen*self.LProp.PopulationNumber*self.LProp.GridBound+pop*self.LProp.GridBound+pos]
        return rate

    # Compute death rate
    def Death(self, pop, gen, pos, N):
        rate = self.DeathPC(pop, gen, pos)
        rate *= N[gen*self.LProp.PopulationNumber*self.LProp.GridBound+pop*self.LProp.GridBound+pos]
        return rate

    # Compute mutation up
    def MuteUp(self, pop, gen, pos, N):
        rate = self.MuteUpPC(pop, gen, pos)
        gen_donor = gen + 1
        while gen_donor<self.LProp.GenBound:
            rate += self.MuteHGTPC(pop, gen, gen_donor, pos, N)
            gen_donor += 1
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute mutation down
    def MuteDown(self, pop, gen, pos, N):
        rate = self.MuteDownPC(pop, gen, pos)
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute move right
    def MoveRight(self, pop, gen, pos, N):
        rate = self.MoveRightPC(pop, gen, pos, N)
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute move left
    def MoveLeft(self, pop, gen, pos, N):
        rate = self.MoveLeftPC(pop, gen, pos, N)
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute switch up
    def SwitchUp(self, pop, gen, pos, N):
        rate = self.SwitchUpPC(pop, gen, pos)
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Compute switch down
    def SwitchDown(self, pop, gen, pos, N):
        rate = self.SwitchDownPC(pop, gen, pos)
        rate *= N[gen * self.LProp.PopulationNumber * self.LProp.GridBound + pop * self.LProp.GridBound + pos]
        return rate

    # Stable wildtype system of equations, N is the ansatz solution vector [N1, N2, N3, ..., Nn, N1, N2, ...]
    def Force(self, N):
        force = np.zeros(self.LProp.GenBound*self.LProp.PopulationNumber*self.LProp.GridBound)
        for gen in range(self.LProp.GenBound):
            for pop in range(self.LProp.PopulationNumber):
                for pos in range(self.LProp.GridBound):
                    live_force = self.Birth(pop, gen, pos, N)-self.Death(pop, gen, pos, N)
                    move_force = 0
                    if pos > 0:
                        move_force += self.MoveRight(pop, gen, pos-1, N)-self.MoveLeft(pop, gen, pos, N)
                    if pos < self.LProp.GridBound-1:
                        move_force += self.MoveLeft(pop, gen, pos+1, N)-self.MoveRight(pop, gen, pos, N)
                    switch_force = 0
                    if pop > 0:
                        switch_force += self.SwitchUp(pop-1, gen, pos, N)-self.SwitchDown(pop, gen, pos, N)
                    if pop<self.LProp.PopulationNumber-1:
                        switch_force += self.SwitchDown(pop+1, gen, pos, N)-self.SwitchUp(pop, gen, pos, N)
                    mut_force=0
                    if gen>0:
                        mut_force += self.MuteUp(pop, gen-1, pos, N)-self.MuteDown(pop, gen, pos, N)
                    if gen<self.LProp.GenBound-1:
                        mut_force += self.MuteDown(pop, gen+1, pos, N)-self.MuteUp(pop, gen, pos, N)
                    force[gen*self.LProp.PopulationNumber*self.LProp.GridBound+pop*self.LProp.GridBound+pos] = live_force+move_force+switch_force+mut_force
        return force

    # Analytical curves
    def Her_adap_rate(self):
        N = self.LProp.CarryingCapacity * (1 - (self.PProp[0].DeathRate + 2 * self.PProp[0].Move) / self.PProp[0].BirthRate)
        kappa = self.PProp[0].MuteUp * N / (2 * (self.PProp[0].DeathRate + 2 * self.PProp[0].Move))
        a = self.PProp[0].MuteUp * N / (self.PProp[0].DeathRate + 2 * self.PProp[0].Move - self.PProp[0].MuteUp)
        c = np.sqrt((self.PProp[0].MuteDown + 2 * self.PProp[0].MuteUp) ** 2 / 4 + self.PProp[0].Move * (self.PProp[0].DeathRate + 2 * self.PProp[0].Move - self.PProp[0].MuteUp))
        if kappa < 1e2:
            wait_time = np.sqrt(np.pi) / (2 * c) * sp.gamma(kappa) / sp.gamma(kappa + 0.5)
        else:
            wait_time = np.sqrt(np.pi / (2 * self.PProp[0].Move * self.PProp[0].MuteUp * N))  # np.sqrt(2*np.pi/a)/(2*c)
        return 1 / wait_time

    # Analytical approximation
    # Compute the cumulative distribution functions F for various paths
    def SD(self, pop, time):
        if self.cD[pop] * time < 1e2:
            sd = (self.cD[pop]*np.exp(self.bD[pop]*time)/(self.cD[pop]*np.cosh(self.cD[pop]*time)+self.bD[pop]*np.sinh(self.cD[pop]*time)))**self.aD[pop]
        else:
            sd = (2*self.cD[pop]/(self.cD[pop]+self.bD[pop]))**self.aD[pop]*np.exp((self.bD[pop]-self.cD[pop])*self.aD[pop]*time)
        return sd

    def FD(self, pop, time):
        n_average = self.MuteUpPC(pop, 0, self.LProp.LD-1)*self.nDeff[pop]*np.tanh(self.cD[pop]*time)/(self.cD[pop]+self.bD[pop]*np.tanh(self.cD[pop]*time))
        fd = self.MoveRightPC(pop, 1, self.LProp.LD-1, self.N)*n_average*self.SD(pop, time)
        return fd

    def SO(self, pop, time):
        so = np.exp(-self.lambdaO[pop]*time)
        return so

    def FO(self, pop, time):
        fo = self.lambdaO[pop]*np.exp(-self.lambdaO[pop]*time)
        return fo

    def SN(self, pop, time):
        if self.cN[pop]*time<1e2:
            sn = (self.cN[pop]*np.exp(self.bN[pop]*time)/(self.cN[pop]*np.cosh(self.cN[pop]*time)+self.bN[pop]*np.sinh(self.cN[pop]*time)))**self.aN[pop]
        else:
            sn = (2*self.cN[pop]/(self.cN[pop]+self.bN[pop]))**self.aN[pop]*np.exp((self.bN[pop]-self.cN[pop])*self.aN[pop]*time)
        return sn

    def FN(self, pop, time):
        n_average = self.MuteUpPC(pop,0, self.LProp.LD+1)*self.nNeff[pop]*np.tanh(self.cN[pop]*time)/(self.cN[pop]+self.bN[pop]*np.tanh(self.cN[pop]*time))
        fn = self.MoveLeftPC(pop,1,self.LProp.LD+1,self.N)*n_average*self.SN(pop,time)
        return fn

    def S(self, time):
        s=1
        for pop in range(self.LProp.PopulationNumber):
            s *= self.SD(pop,time) * self.SO(pop,time) * self.SN(pop,time)
        return s

    def F(self, time):
        f=0
        for pop in range(self.LProp.PopulationNumber):
            termD = 1
            termO = 1
            termN = 1
            for popl in range(self.LProp.PopulationNumber):
                if popl == pop:
                    termD *= self.FD(popl,time) * self.SO(popl,time) * self.SN(popl,time)
                    termO *= self.SD(popl, time) * self.FO(popl, time) * self.SN(popl, time)
                    termN *= self.SD(popl, time) * self.SO(popl, time) * self.FN(popl, time)
                else:
                    termD *= self.SD(popl, time) * self.SO(popl, time) * self.SN(popl, time)
                    termO *= self.SD(popl, time) * self.SO(popl, time) * self.SN(popl, time)
                    termN *= self.SD(popl, time) * self.SO(popl, time) * self.SN(popl, time)
            f += termD + termO + termN
        return f

    # Find the adaptation rate by integrating
    def integrand(self, time):
        return time * self.F(time)

    # Computes evolutionary time
    def T1(self):
        if self.N[self.LProp.GridBound-1] > 0:
            t = integrate.quad(self.integrand, 0, np.inf)
            t1 = t[0]
        else:
            t1 = np.inf
        return t1

    # Integration by Runge-Kutta
    def RK4_step(self, h, x):
        k1 = self.Force(x) * h
        k2 = self.Force(x + k1 / 2) * h
        k3 = self.Force(x + k2 / 2) * h
        k4 = self.Force(x + k3) * h
        x_new = x + (k1 + k2 * 2 + k3 * 2 + k4) / 6
        return x_new

    # Checks if some wild-type is outcompeted for some mutants in the overlap regions
    def not_outcompeted(self, x):
        not_out = True
        for pop in range(self.LProp.PopulationNumber):
            if x[pop*self.LProp.GridBound+self.LProp.LD]<x[self.LProp.GridBound*self.LProp.PopulationNumber+pop*self.LProp.GridBound+self.LProp.LD]:
                    not_out = False
        return not_out

    # Checks if overall wild-type is outcompeted by overall mutants for some population
    def overall_not_outcompeted(self, x):
        not_out = True
        for pop in range(self.LProp.PopulationNumber):
            wild_type=0
            mutants=0
            for pos in range(self.LProp.GridBound):
                wild_type += x[pop*self.LProp.GridBound+pos]
                mutants += x[self.LProp.GridBound*self.LProp.PopulationNumber+pop*self.LProp.GridBound+pos]
            if wild_type<mutants:
                    not_out = False
        return not_out

    # Computes ecological time
    def D(self, h, t_max):
        x = np.zeros(2*self.LProp.GridBound*self.LProp.PopulationNumber)
        for pop in range(self.LProp.PopulationNumber):
            for pos in range(self.LProp.GridBound):
                if pos < self.LProp.LD:
                    x[pop*self.LProp.GridBound+pos] = self.N[pop*self.LProp.GridBound+pos]-self.MuteUpPC(pop,0,pos)*self.N[pop*self.LProp.GridBound+pos]/(self.bD[pop]+self.cD[pop])
                    x[self.LProp.GridBound*self.LProp.PopulationNumber+pop*self.LProp.GridBound+pos] = self.MuteUpPC(pop,0,pos)*self.N[pop*self.LProp.GridBound+pos]/(self.bD[pop]+self.cD[pop])
                elif pos == self.LProp.LD:
                    x[pop*self.LProp.GridBound+pos] = self.N[pop*self.LProp.GridBound+pos] - 1
                    x[self.LProp.GridBound*self.LProp.PopulationNumber+pop*self.LProp.GridBound+pos] = 1
                else:
                    x[pop*self.LProp.GridBound+pos] = self.N[pop*self.LProp.GridBound+pos]-self.MuteUpPC(pop,0,pos)*self.N[pop*self.LProp.GridBound+pos]/(self.bN[pop]+self.cN[pop])
                    x[self.LProp.GridBound*self.LProp.PopulationNumber+pop*self.LProp.GridBound+pos] = self.MuteUpPC(pop,0,pos)*self.N[pop*self.LProp.GridBound+pos]/(self.bN[pop]+self.cN[pop])
        D = 0
        while self.not_outcompeted(x) and (D < t_max):
            x = self.RK4_step(h, x)
            D += h
        return D

    # computes usual adaptation rate
    def adap_rate(self, h, t_max):
        T1 = self.T1()
        D = self.D(h, t_max)
        print('T1=' + str(T1))
        print('D=' + str(D))
        if np.isinf(T1) or T1 < 0 or T1+D == 0:
            return 'NaN'
        else:
            return 1/(T1+D)

    # computes the adaptation rate where D is computed from LD-1, this appears in simulations
    def retarded_adap_rate(self, h, t_max):
        T1 = self.T1()
        if self.LProp.LD > 1:
            self.LProp.LD -= 1
            self.Initialize()
            D = self.D(h, t_max)
        else:
            D = 0
        print('T1=' + str(T1))
        print('D=' + str(D))
        if np.isinf(T1) or T1 < 0 or T1 + D == 0:
            return 'NaN'
        else:
            return 1 / (T1 + D)

    # computes overall deterministic waiting time starting from x_init (2*GridBound*PopNum entries), time step h, up to t_max
    def deterministic_time(self, h, t_max, x_init):
        x = x_init
        D = 0
        while self.overall_not_outcompeted(x) and (D < t_max):
            x = self.RK4_step(h, x)
            D += h
        return D

    # computes swarming tolerance of swarm_pop
    def swarming_tolerance(self, swarm_pop):
        # find the maximal position which allows swarming
        max_pos = 0
        for pos in range(self.LProp.GridBound):
            if self.N[swarm_pop * self.LProp.GridBound + pos] >= self.PProp[swarm_pop].SwarmingDensity:
                max_pos = pos
        # find swarming tolerance
        return max_pos+1-self.LProp.LD

    # computes swarming phenotypic waiting times
    # INPUT: innoc_N = [for each pop] - located in the first compartment; swarm_pop - measured swarming population
    # OUTPUT: wait_times = [times when reached swarming at allowed positions]
    def phen_wait_time(self, h, t_max, innoc_N, swarm_pop):
        # neglect mutations
        self.LProp.GenBound = 1
        # set initial conditions
        x = np.zeros(self.LProp.GridBound*self.LProp.PopulationNumber)
        for pop in range(self.LProp.PopulationNumber):
            x[pop * self.LProp.GridBound] = innoc_N[pop]
        time = 0
        # find the maximal position which allows swarming
        max_pos = 0
        for pos in range(self.LProp.GridBound):
            if self.N[swarm_pop*self.LProp.GridBound+pos]>=self.PProp[swarm_pop].SwarmingDensity:
                max_pos = pos
        wait_times = np.ones(max_pos+1)*np.inf
        # do the simulation
        adap_pos = 0
        while adap_pos <= max_pos and (time < t_max):
            x = self.RK4_step(h, x)
            time += h
            if x[swarm_pop * self.LProp.GridBound+adap_pos]>=self.PProp[swarm_pop].SwarmingDensity:
                wait_times[adap_pos] = time
                adap_pos += 1
        # put mutations back
        self.LProp.GenBound = 2
        # return the result
        return wait_times

    # computes swarming phenotypic adaptation rate from position pos to pos+1, type=just swarm(0), need to reproduce(1)
    def phen_adap_rate(self, pos, h, t_max, innoc_N, swarm_pop):
        wait_times = self.phen_wait_time(h, t_max, innoc_N, swarm_pop)
        if pos < len(wait_times)-1:
            if np.isinf(wait_times[pos+1]):
                print("Warning! Prolong simulation time to find when swarming occurs at "+str(pos+1))
                adap_rate = 0
            elif wait_times[pos+1] == wait_times[pos]:
                print("Warning! The adaptation seems to be instanteneous. Decrease simulation time-step.")
                adap_rate = np.inf
            else:
                adap_rate = 1/(wait_times[pos+1]-wait_times[pos])
        else:
            print("Warning! Need genotypic mutation to swarm at "+str(pos+1))
            adap_rate = 0
        return adap_rate

    # computes adaptation rate and treatment efficacy, outputs [adap_rate, efficacy]
    def adap_efficacy(self, h, t_max):
        # find population with maximal net birth rate (this would naturally outcompete the others)
        pop_max = 0
        net_birth = 0
        for pop in range(self.LProp.PopulationNumber):
            new_net_birth = self.PProp[pop].BirthRate - self.PProp[pop].DeathRate
            if new_net_birth > net_birth:
                pop_max = pop
                net_birth = new_net_birth
        # compute size of population in absence of antibiotics
        N0 = self.LProp.GridBound*self.LProp.CarryingCapacity*(1-self.PProp[pop_max].DeathRate/self.PProp[pop_max].BirthRate)
        # compute adaptation rate
        adap_rate = self.adap_rate(h, t_max)
        # compute total stable wild-type and treatment efficacy
        if self.N[0]>0:
            N_tot = sum(self.N)
            eff = (N0 - N_tot) / N0
            if not adap_rate=='NaN':
                eff *= np.tanh(self.PProp[pop_max].DeathRate / adap_rate)
        else:
            eff = 1
        # return efficacy
        return [adap_rate, eff]

# Find Critical Motility
def CritMot(nu_max, error, pop, p_prop, l_prop):
    p_prop[pop].Move = nu_max
    stair = Staircase(l_prop, p_prop)
    if stair.N[0] > 0:
        return None
    else:
        nu_min = p_prop[pop].DeathRate
        while nu_max - nu_min > error:
            nu = (nu_min + nu_max) / 2
            p_prop[pop].Move = nu
            stair = Staircase(l_prop, p_prop)
            if stair.N[0] > 0:
                nu_min = nu
            else:
                nu_max = nu
        return nu