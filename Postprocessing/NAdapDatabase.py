###
# IMPORTS
###
import numpy as np
from os import path
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
import copy
import NStaircase as sc
from io import TextIOWrapper
import zipfile


###
# MAIN CLASS
###

# Experiment List (A=oldest version, B=newest version):
# 0=MotSwitchA -> params=(move0, move1, switch)
# 1=ChemotaxA -> params=(move, chemotax)
# 2=StandardA -> params=(move, death)
# 3=SpecialA -> params=(move, type), type=0(resist cost),=1(sbirth),=2(sdeath)
# 4=ResistCostA -> params=(move, resist_cost)
# 5=StressBirthA -> params=(move, stress_birth)
# 6=StressDeathA -> params=(move, stress_death)
# 7=SpecialB -> params=(move, type), type=0(resist cost),=1(sbirth),=2(sbirth,sdeath)
# 8=ResistCostB -> params=(move, resist_cost)
# 9=StressBirthDeathA -> params=(move, stress_death)
# 10=StressDeathB -> params=(move, stress_death)
# 11=SwitchA -> params=(switch, death)
# 12=MotSwitchDeadlyA -> params=(move0, move1, switch)
# 13=SwarmLagA -> params=(swarm_move, init_cell)
# 14=SwarmDeadlyA -> params=(swarm_move, swarm_density, death)
# 15=HGTA -> params=(move, hgt)
# 16=HGTDeadlyA -> params=(move, hgt)
# 17=DensSwitchA -> params=(death, switch_density, move_low_dens, move_high_dens)
# 18=ResCostA -> params=(move, resist_cost)
# 19=ResCostDeadlyA -> params=(move, resist_cost)

class AdapDatabase:

    ###
    # DATA ACCESS
    ###

    # input: folder_paths = ['Fig2/FigIIv1/', ...], types of the folder_paths = [0 = MotSwitch, 1=Chemotax, 2=Standard]
    #        zip = None (usual folders), list of zip file names per each folder path
    def __init__(self, folder_paths, types, zip):
        self.folder_paths = folder_paths
        if zip is None:
            self.zip = False
        else:
            self.zip = True
            self.zip_files = []
            for zip_name in zip:
                self.zip_files.append(zipfile.ZipFile(zip_name))
        possible_experiment_types = [MotSwitchA(), ChemotaxA(), StandardA(), SpecialA(), ResistCostA(), StressBirthA(), StressDeathA(), SpecialB(), ResistCostB(), StressBirthDeathA(), StressDeathB(), SwitchA(), MotSwitchDeadlyA(), SwarmLagA(), SwarmDeadlyA(), HGTA(), HGTDeadlyA(), DensSwitchA(), ResCostA(), ResCostDeadlyA()]
        self.exp_type = []
        for type in types:
            self.exp_type.append(copy.deepcopy(possible_experiment_types[type]))

    # read snapshot file
    # output: [Time,GenSpaceTot[pop][gen][pos]]
    # input: sim_id(=0,...), task_id, population evolving mutation, LD, time_type (='T','D','PA') -> 'OutAdap_taskID_ID_pop_LDT/D.csv'
    def read_gen_space_tot(self, sim_id, task_id, population, LD, time_type):
        time_types = ['T', 'D', 'PA']
        file_stem = "OutAdap_"+str(int(task_id))+"_"+str(self.exp_type[sim_id].taskIDtoID(task_id))
        snapshot_file = self.folder_paths[sim_id] + file_stem + "_"+str(population)+"_"+str(LD)+time_types[time_type]+".csv"
        GenSpaceTot = np.zeros((self.exp_type[sim_id].PopulationNumber, self.exp_type[sim_id].GenBound, self.exp_type[sim_id].GridBound), dtype=float)
        if (not self.zip and path.exists(snapshot_file)) or (self.zip and (snapshot_file in self.zip_files[sim_id].namelist())):
            if not self.zip:
                file = open(snapshot_file, 'r')
            else:
                file = TextIOWrapper(self.zip_files[sim_id].open(snapshot_file, 'r'), 'utf-8')
            reader = list(csv.reader(file, skipinitialspace=True))
            if len(reader) != 0:
                Time = float(reader[1][1])
                for pop in range(self.exp_type[sim_id].PopulationNumber):
                    for gen in range(self.exp_type[sim_id].GenBound):
                        for pos in range(self.exp_type[sim_id].GridBound):
                            GenSpaceTot[pop][gen][pos] = float(reader[pop*(self.exp_type[sim_id].GenBound+1)+gen+3][pos])
                return [Time, GenSpaceTot]
            else:
                print("Warning! Tried to read gen_space_tot but the snapshot_file is empty: " + snapshot_file)
                return None
        else:
            print("Warning! Tried to read gen_space_tot but the snapshot_file does not exist: " + snapshot_file)
            return None

    # read times file
    # output: [TTimes[population][LD], DTimes[population][LD], PATimes[population][LD]]
    # notice: LD=0,...,GridBound-1 with TTimes/DTimes waiting times to adapt in LD+1, PATimes to swarm in LD; unused entries filled with 0s
    # input: sim_id, task_id
    def read_times(self, sim_id, task_id):
        file_stem = "OutAdap_" + str(int(task_id)) + "_" + str(self.exp_type[sim_id].taskIDtoID(task_id))
        times_file = self.folder_paths[sim_id] + file_stem + "_Times.csv"
        if (not self.zip and path.exists(times_file)) or (self.zip and (times_file in self.zip_files[sim_id].namelist())):
            if not self.zip:
                file = open(times_file, 'r')
            else:
                file = TextIOWrapper(self.zip_files[sim_id].open(times_file, 'r'), 'utf-8')
            reader = list(csv.reader(file, skipinitialspace=True))
            if len(reader) != 0:
                grid_bound = self.exp_type[sim_id].GridBound
                TTimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
                DTimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
                PATimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
                num_read_lines = 0
                for pop in range(self.exp_type[sim_id].PopulationNumber):
                    for LD in range(grid_bound):
                        if LD < len(reader[num_read_lines+2]) and reader[num_read_lines+2][LD] != '':
                            TTimes[pop][LD] = float(reader[num_read_lines+2][LD])
                        if LD < len(reader[num_read_lines+4]) and reader[num_read_lines+2][LD] != '':
                            DTimes[pop][LD] = float(reader[num_read_lines+4][LD])
                        if self.exp_type[sim_id].ConsiderSwarm[pop]:
                            if LD < len(reader[num_read_lines + 6]) and reader[num_read_lines+2][LD] != '':
                                PATimes[pop][LD] = float(reader[num_read_lines + 6][LD])
                    if not self.exp_type[sim_id].ConsiderSwarm[pop]:
                        num_read_lines += 5
                    else:
                        num_read_lines += 7
                return [TTimes, DTimes, PATimes]
            else:
                print("Warning! Tried to read times but the times_file is empty: " + times_file)
                return None
        else:
            print("Warning! Tried to read times but the times_file does not exist: " + times_file)
            print("Instead, reading from GenSpaceTotFiles will be attempted.")
            grid_bound = self.exp_type[sim_id].GridBound
            adap_num = self.exp_type[sim_id].adap_num
            TTimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
            DTimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
            PATimes = np.zeros((self.exp_type[sim_id].PopulationNumber, grid_bound))
            output = False
            for pop in range(self.exp_type[sim_id].PopulationNumber):
                for LD in range(adap_num):
                    if not self.exp_type[sim_id].ConsiderSwarm[pop]:
                        result = self.read_gen_space_tot(sim_id, task_id, pop, LD+1, 0)
                        if result is not None:
                            output = True
                            TTimes[pop][LD] = result[0]
                        result = self.read_gen_space_tot(sim_id, task_id, pop, LD+1, 1)
                        if result is not None:
                            output = True
                            DTimes[pop][LD] = result[0]
                if self.exp_type[sim_id].ConsiderSwarm[pop]:
                    for LD in range(grid_bound):
                        result = self.read_gen_space_tot(sim_id, task_id, pop, LD, 2)
                        if result is not None:
                            output = True
                            PATimes[pop][LD] = result[0]
            if output is True:
                return [TTimes, DTimes, PATimes]
            else:
                return None

    # compute average adaptation rate (all experiments must have the same GridBound)
    # output: [adap[type][LD] (type, 0=T, 1=D, 2=PA), standard_deviation[type][LD]]
    # notice: LD=0,...,GridBound-1 with T/D to adapt in LD+1, PA to swarm in LD; unused entries filled with 0s
    # use harmonic means + standard deviation as in:
    # https://stats.stackexchange.com/questions/7471/can-the-standard-deviation-be-calculated-for-harmonic-mean
    # input: params (either a single tuple - used for all folders, or a list of tuples for each folder)
    def adap_rate(self, params):
        # change a single tuple to a list of tuples
        if not isinstance(params, list):
            params_list = [params for sim_id in range(len(self.exp_type))]
            params = params_list
        # find the total number of repeats, and check uniformity of GridBound
        grid_bound = self.exp_type[0].GridBound
        total_repeats = 0
        for sim_id in range(len(self.exp_type)):
            total_repeats += self.exp_type[sim_id].repeat
            if grid_bound != self.exp_type[sim_id].GridBound:
                print("GridBounds in adap_rate computation differ. Function adap_rate might misbehave!")
        # prepare containers
        adap = np.zeros((3, grid_bound))                      # adaptation rate
        standard_deviation = np.zeros((3, grid_bound))        # standard deviation
        delta = np.zeros((3, grid_bound, total_repeats))      # time differences between steps, sum over all computation
        num = np.zeros((3, grid_bound), dtype=int)            # number of time differences added to the sum
        for sim_id in range(len(self.exp_type)):
            if params[sim_id] is not None:
                for run_index in range(self.exp_type[sim_id].repeat):
                    taskID = self.exp_type[sim_id].IDtoTaskID(params[sim_id], run_index)
                    times = self.read_times(sim_id, taskID)     # first index 0=T, 1=D, 2=PA
                    if times is not None:
                        for type in range(3):
                            lower_t = 0
                            for LD in range(grid_bound):
                                upper_t = times[type][0][LD]  # find minimal non-zero element over all populations
                                for pop in range(len(times[type])):
                                    if upper_t*times[type][pop][LD] == 0:  # if one of them is zero
                                        upper_t = upper_t+times[type][pop][LD]  # choose the non-zero
                                    else:
                                        upper_t = min(upper_t, times[type][pop][LD])    # else choose minimum
                                if upper_t > 0:
                                    if (LD == 0) or (lower_t > 0):
                                        delta[type][LD][num[type][LD]] = upper_t-lower_t
                                        num[type][LD] += 1
                                lower_t = upper_t
        for type in range(3):
            for LD in range(grid_bound):
                if num[type][LD] > 0:
                    expected_time = np.sum(delta[type][LD][:num[type][LD]])/num[type][LD]
                    time_variance = np.sum((delta[type][LD][:num[type][LD]]-expected_time)**2)/num[type][LD]
                    if expected_time > 0:
                        standard_deviation[type][LD] = np.sqrt(time_variance/(num[type][LD]*expected_time**4))
                        adap[type][LD] = 1/expected_time
                    else:
                        standard_deviation[type][LD] = 'NaN'
                        adap[type][LD] = np.inf
                else:
                    standard_deviation[type][LD] = 'NaN'
                    adap[type][LD] = 'NaN'
        return [adap, standard_deviation]

    # compute swarming tolerance, expected number of compartments under the staircase to swarm in before mutants outcompete wildtype gen
    # notice: all experiments must have the same GridBound, GenBound!
    # output: [swarm_tol[gen], standard_deviation[gen]], gen=0,1,...,GenBound-1
    # input: params (either a single tuple - used for all folders, or a list of tuples for each folder), pop (swarming population)
    def swarm_tolerance(self, params, pop):
        # change a single tuple to a list of tuples
        if not isinstance(params, list):
            params_list = [params for sim_id in range(len(self.exp_type))]
            params = params_list
        # find the total number of repeats, and check uniformity of GridBound
        grid_bound = self.exp_type[0].GridBound
        gen_bound = self.exp_type[0].GenBound
        total_repeats = 0
        for sim_id in range(len(self.exp_type)):
            total_repeats += self.exp_type[sim_id].repeat
            if grid_bound != self.exp_type[sim_id].GridBound:
                print("GridBounds in adap_rate computation differ. Function adap_rate might misbehave!")
            if gen_bound != self.exp_type[sim_id].GenBound:
                print("GenBounds in adap_rate computation differ. Function adap_rate might misbehave!")
        # prepare containers
        swarm_tol = np.zeros(gen_bound)                   # swarming tolerance
        standard_deviation = np.zeros(gen_bound)         # standard deviation
        swarm_x = np.zeros((gen_bound, total_repeats))      # swarming tolerances in particular experiments
        num = np.zeros(gen_bound, dtype=int)            # number of time differences added to the sum
        for sim_id in range(len(self.exp_type)):
            if params[sim_id] is not None:
                for run_index in range(self.exp_type[sim_id].repeat):
                    taskID = self.exp_type[sim_id].IDtoTaskID(params[sim_id], run_index)
                    times = self.read_times(sim_id, taskID)     # first index 0=T, 1=D, 2=PA
                    if times is not None:
                        for gen in range(gen_bound):
                            D = times[1][pop][gen]
                            if D>0:
                                swarm_pos = 0
                                while (swarm_pos<grid_bound-gen-1) and (times[2][pop][swarm_pos+1+gen]<D) and (times[2][pop][swarm_pos+1+gen] != 0):
                                    swarm_pos += 1
                                swarm_x[gen][num[gen]] = swarm_pos
                                num[gen] += 1
        for gen in range(gen_bound):
            if num[gen] > 1:
                swarm_tol[gen] = np.sum(swarm_x[gen][:num[gen]])/num[gen]
                standard_deviation[gen] = np.sqrt(np.sum((swarm_x[gen][:num[gen]]-swarm_tol[gen])**2)/(num[gen]-1))
            else:
                swarm_tol[gen] = 'NaN'
                standard_deviation[gen] = 'NaN'
        return [swarm_tol, standard_deviation]

    # compute average waiting time (all experiments must have the same GridBound)
    # output: [wait_time[type][LD] (type, 0=T, 1=D, 2=PA), standard_deviation[type][LD]]
    # notice: LD=0,...,GridBound-1 with T/D to adapt in LD+1, PA to swarm in LD; unused entries filled with 0s
    # input: params (either a single tuple - used for all folders, or a list of tuples for each folder)
    def waiting_time(self, params):
        # change a single tuple to a list of tuples
        if not isinstance(params, list):
            params_list = [params for sim_id in range(len(self.exp_type))]
            params = params_list
        # find the total number of repeats, and check uniformity of GridBound
        grid_bound = self.exp_type[0].GridBound
        total_repeats = 0
        for sim_id in range(len(self.exp_type)):
            total_repeats += self.exp_type[sim_id].repeat
            if grid_bound != self.exp_type[sim_id].GridBound:
                print("GridBounds in adap_rate computation differ. Function adap_rate might misbehave!")
        # prepare containers
        wait_time = np.zeros((3, grid_bound))                      # adaptation rate
        standard_deviation = np.zeros((3, grid_bound))        # standard deviation
        delta = np.zeros((3, grid_bound, total_repeats))      # time differences between steps, sum over all computation
        num = np.zeros((3, grid_bound), dtype=int)            # number of time differences added to the sum
        for sim_id in range(len(self.exp_type)):
            if params[sim_id] is not None:
                for run_index in range(self.exp_type[sim_id].repeat):
                    taskID = self.exp_type[sim_id].IDtoTaskID(params[sim_id], run_index)
                    times = self.read_times(sim_id, taskID)     # first index 0=T, 1=D, 2=PA
                    if times is not None:
                        for type in range(3):
                            lower_t = 0
                            for LD in range(grid_bound):
                                upper_t = times[type][0][LD]  # find minimal non-zero element over all populations
                                for pop in range(len(times[type])):
                                    if upper_t*times[type][pop][LD] == 0:  # if one of them is zero
                                        upper_t = upper_t+times[type][pop][LD]  # choose the non-zero
                                    else:
                                        upper_t = min(upper_t, times[type][pop][LD])    # else choose minimum
                                if upper_t > 0:
                                    if (LD == 0) or (lower_t > 0):
                                        delta[type][LD][num[type][LD]] = upper_t-lower_t
                                        num[type][LD] += 1
                                lower_t = upper_t
        for type in range(3):
            for LD in range(grid_bound):
                if num[type][LD] > 1:
                    wait_time[type][LD] = np.sum(delta[type][LD][:num[type][LD]])/num[type][LD]
                    standard_deviation[type][LD] = np.sqrt(np.sum((delta[type][LD][:num[type][LD]]-wait_time[type][LD])**2)/(num[type][LD]-1))
                else:
                    standard_deviation[type][LD] = 'NaN'
                    wait_time[type][LD] = 'NaN'
        return [wait_time, standard_deviation]

    ###
    # PLOTTING
    ###

    # make the motility-adaptation rate plot of given time_type, LD, and death parameter
    # output: [min_adap_rate, max_adap_rate, plots]
    # input: death (must match existing death_param), time_type (0=T, 1=D), LD,
    #        analytics (t/f, poin_num, precision, Hermsen, hom_sel), ax, min_max, labels (t/f), legend_name (string), background (t/f)
    # note: assume all input data are of the same StandardA experiment
    def standard_motility_adap_plot(self, death, LD, time_type, analytics, ax, min_max, labels, legend_name, background):
        # prepare adaptation rates, and find minimal and maximal adaptation rate
        plots = []
        move_num = len(self.exp_type[0].move_params)
        adap_rates = np.zeros(move_num)
        crit_mot = None
        for move_index in range(move_num):
            move = self.exp_type[0].move_params[move_index]
            # get adaptation rate
            if self.adap_rate((move, death)) is None:
                adap_rate = 'NaN'
            else:
                adap_rate = self.adap_rate((move, death))[0][time_type][LD - 1]
                if np.isnan(adap_rate):
                    adap_rate = 'NaN'
            adap_rates[move_index] = adap_rate
            # get critical motility
            if crit_mot is None:
                if adap_rate == 'NaN':
                    crit_mot = move
            else:
                if adap_rate != 'NaN':
                    crit_mot = None
        if min_max is None:
            max = adap_rates.max()
            min = np.nanmin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
        # make the plot
        p = ax.scatter(self.exp_type[0].move_params, adap_rates, label=legend_name[0])
        plots.append(p)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(self.exp_type[0].move_params[0] / 10, self.exp_type[0].move_params[-1] * 10)
        # make the background
        if background:
            if crit_mot is None:
                ax.axvspan(self.exp_type[0].move_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, self.exp_type[0].move_params[-1] * 10, facecolor="tab:orange", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.5*min, "$\\delta$", fontsize=12, color='black')
                ax.text(self.exp_type[0].move_params[0] / 5, 4 * min, "Low Motility", fontsize=15, color='green', rotation=90)
                ax.text(10 * death, 4 * min, "High Motility", fontsize=15, color='orange', rotation=90)
            else:
                ax.axvspan(self.exp_type[0].move_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, crit_mot, facecolor="tab:orange", alpha=0.3)
                ax.axvspan(crit_mot, self.exp_type[0].move_params[-1] * 10, facecolor="tab:red", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.5*min, "$\\delta$", fontsize=12, color='black')
                ax.axvline(crit_mot, color='black', linestyle="--")
                ax.text(1.2 * crit_mot, 1.5*min, "$\\nu_c$", fontsize=12)
                ax.text(self.exp_type[0].move_params[0] / 5, 4 * min, "Low Motility", fontsize=15, color='green', rotation=90)
                ax.text(2 * death, 4 * min, "High Motility", fontsize=15, color='orange', rotation=90)
                ax.text(2 * crit_mot, 4 * min, "Deadly Motility", fontsize=15, color='red', rotation=90)
        # plot the labels
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)")
            ax.set_ylabel("adaptation $a_R$ (1/h)")
        # plot the analytics
        if analytics[0]:
            point_num = analytics[1]
            move_params = np.logspace(np.log10(self.exp_type[0].move_params[0]/2), np.log10(self.exp_type[0].move_params[-1]*2), num=point_num)
            adap_params = np.zeros(point_num)
            if len(analytics)>4 and analytics[4]:
                ax_twin = ax.twinx()
                sel_params = np.zeros(point_num)
                ax_twin.set_yscale("log")
            if analytics[3]:
                hermsen_params = np.zeros(point_num)
            l_prop = sc.LatProp()
            l_prop.LD = LD
            p_prop = [sc.PopProp()]
            p_prop[0].DeathRate = death
            for i in range(point_num):
                p_prop[0].Move = move_params[i]
                stair = sc.Staircase(l_prop, p_prop)
                if analytics[3]:
                    hermsen_params[i] = stair.Her_adap_rate()
                adap_params[i] = stair.retarded_adap_rate(analytics[2], 1e4)
                if len(analytics) > 4 and analytics[4]:
                    sel_params[i] = stair.homog_selection(analytics[2], 1e4)
            if crit_mot is None:
                mask = np.isfinite(adap_params)
                p = ax.plot(move_params[mask], adap_params[mask], label=legend_name[1])
                if len(analytics) > 4 and analytics[4]:
                    mask = np.isfinite(sel_params)
                    q = ax_twin.plot(move_params[mask], sel_params[mask], label=legend_name[1], color='red')
            else:
                mask1 = np.zeros(point_num, dtype=bool)
                mask2 = np.zeros(point_num, dtype=bool)
                for i in range(point_num):
                    if np.isfinite(adap_params[i]):
                        if move_params[i] < crit_mot:
                            mask1[i] = True
                        else:
                            mask2[i] = True
                p = ax.plot(move_params[mask1], adap_params[mask1], label=legend_name[1])
                ax.plot(move_params[mask2], adap_params[mask2], '--', color=p[0].get_color(), label=legend_name[1])
                if len(analytics) > 4 and analytics[4]:
                    mask1 = np.zeros(point_num, dtype=bool)
                    mask2 = np.zeros(point_num, dtype=bool)
                    for i in range(point_num):
                        if np.isfinite(sel_params[i]):
                            if move_params[i] < crit_mot:
                                mask1[i] = True
                            else:
                                mask2[i] = True
                    q = ax_twin.plot(move_params[mask1], sel_params[mask1], label=legend_name[1], color='red')
                    ax.plot(move_params[mask2], sel_params[mask2], '--', color=p[0].get_color(), label=legend_name[1])
            plots.append(p[0])
            if len(analytics) > 4 and analytics[4]:
                plots.append(q[0])
                ax_twin.set_ylabel("selection $s$", fontsize=12)
                #ax_twin.set_ylim(2e-3, 30)
            if analytics[3]:
                p = ax.plot(move_params, hermsen_params, '--', color='grey', label="Hermsen, 2012")
                plots.append(p[0])
        ax.set_ylim(min, max)
        return [min, max, plots]

    # make the switching-adaptation rate plot of given time_type, LD, and death parameter
    # output: [min_adap_rate, max_adap_rate]
    # input: death (must match existing death_param), time_type (0=T, 1=D), LD, asymptotes (t/f)
    #        analytics (t/f, poin_num, precision), ax, min_max, labels (t/f), legend_name (string), background (t/f)
    # note: assume all input data are of the same SwitchA experiment
    def switch_adap_plot(self, death, LD, time_type, analytics, asymptotes, ax, min_max, labels, legend_name, background):
        # prepare adaptation rates, and find minimal and maximal adaptation rate
        plots = []
        switch_num = len(self.exp_type[0].switch_params)
        adap_rates = np.zeros(switch_num)
        crit_switch = None
        for switch_index in range(switch_num):
            switch = self.exp_type[0].switch_params[switch_index]
            # get adaptation rate
            if self.adap_rate((switch, death)) is None:
                adap_rate = 'NaN'
            else:
                adap_rate = self.adap_rate((switch, death))[0][time_type][LD - 1]
                if np.isnan(adap_rate):
                    adap_rate = 'NaN'
            adap_rates[switch_index] = adap_rate
            # get critical motility
            if crit_switch is None:
                if adap_rate == 'NaN':
                    crit_switch = switch
            else:
                if adap_rate != 'NaN':
                    crit_switch = None
        if min_max is None:
            max = adap_rates.max()
            min = np.nanmin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
        # make the plot
        p = ax.scatter(self.exp_type[0].switch_params, adap_rates, label=legend_name[0])
        plots.append(p)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(self.exp_type[0].switch_params[0] / 10, self.exp_type[0].switch_params[-1] * 10)
        # make the background
        if background:
            if crit_switch is None:
                ax.axvspan(self.exp_type[0].switch_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, self.exp_type[0].switch_params[-1] * 10, facecolor="tab:orange", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.5*min, "$\\delta$", fontsize=12, color='black')
                ax.text(self.exp_type[0].switch_params[0] / 5, 1.5 * min, "Low Switching", fontsize=15, color='green', rotation=90)
                ax.text(2 * death, 1.5 * min, "High Switching", fontsize=15, color='orange', rotation=90)
            else:
                ax.axvspan(self.exp_type[0].switch_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, crit_switch, facecolor="tab:orange", alpha=0.3)
                ax.axvspan(crit_switch, self.exp_type[0].switch_params[-1] * 10, facecolor="tab:red", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.5*min, "$\\delta$", fontsize=12, color='black')
                ax.axvline(crit_switch, color='black', linestyle="--")
                ax.text(1.2 * crit_switch, 1.5*min, "$\\nu_c$", fontsize=12)
                ax.text(self.exp_type[0].switch_params[0] / 5, 1.5 * min, "Low Switching", fontsize=15, color='green', rotation=90)
                ax.text(2 * death, 1.5 * min, "High Switching", fontsize=15, color='orange', rotation=90)
                ax.text(4 * crit_switch, 1.5 * min, "Deadly Switching", fontsize=15, color='red', rotation=90)
        # plot the labels
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)")
            ax.set_ylabel("adaptation $a_2$ (1/h)")
        # plot low switching asymptote
        if asymptotes[0]:
            switch_params = [self.exp_type[0].switch_params[0]/10, death]
            l_prop = sc.LatProp()
            l_prop.LD = LD
            p_prop = [sc.PopProp()]
            p_prop[0].Move = self.exp_type[0].Move[0]
            p_prop[0].DeathRate = death
            stair = sc.Staircase(l_prop, p_prop)
            adap_rate = stair.retarded_adap_rate(1e-2, 1e4)
            adap_params = [adap_rate, adap_rate]
            p = ax.plot(switch_params, adap_params, '--', color='grey')
            plots.append(p[0])
        # plot high switching asymptote
        if asymptotes[1]:
            if crit_switch is None:
                switch_params = [death, self.exp_type[0].switch_params[-1] * 10]
                l_prop = sc.LatProp()
                l_prop.LD = LD
                p_prop = [sc.PopProp()]
                p_prop[0].Move = (self.exp_type[0].Move[1]+self.exp_type[0].Move[0])/2
                p_prop[0].DeathRate = death
                stair = sc.Staircase(l_prop, p_prop)
                adap_rate = stair.retarded_adap_rate(1e-2, 1e4)
                adap_params = [adap_rate, adap_rate]
                p = ax.plot(switch_params, adap_params, '--', color=plots[0].get_facecolor())
                plots.append(p[0])
        # plot the analytics
        if analytics[0]:
            point_num = analytics[1]
            switch_params = np.logspace(np.log10(self.exp_type[0].switch_params[0]/2), np.log10(self.exp_type[0].switch_params[-1]*2), num=point_num)
            adap_params = np.zeros(point_num)
            l_prop = sc.LatProp()
            l_prop.LD = LD
            l_prop.PopulationNumber = 2
            p_prop = [sc.PopProp(), sc.PopProp()]
            p_prop[0].Move = self.exp_type[0].Move[0]
            p_prop[1].Move = self.exp_type[0].Move[1]
            for pop in range(2):
                p_prop[pop].DeathRate = death
            for i in range(point_num):
                for pop in range(2):
                    p_prop[pop].SwitchUp = switch_params[i]
                    p_prop[pop].SwitchDown = switch_params[i]
                stair = sc.Staircase(l_prop, p_prop)
                adap_params[i] = stair.retarded_adap_rate(analytics[2], 1e4)
            if crit_switch is None:
                mask = np.isfinite(adap_params)
                p = ax.plot(switch_params[mask], adap_params[mask], label=legend_name[1])
            else:
                mask1 = np.zeros(point_num, dtype=bool)
                mask2 = np.zeros(point_num, dtype=bool)
                for i in range(point_num):
                    if np.isfinite(adap_params[i]):
                        if switch_params[i] < crit_switch:
                            mask1[i] = True
                        else:
                            mask2[i] = True
                p = ax.plot(switch_params[mask1], adap_params[mask1], label=legend_name[1])
                ax.plot(switch_params[mask2], adap_params[mask2], '--', color=p[0].get_color(), label=legend_name[1])
            plots.append(p[0])
        ax.set_ylim(min, max)
        return [min, max, plots]

    # make the heatmap plot of adaptation rat of given time_type, LD
    # output: [min_adap_rate, max_adap_rate, im]
    # input: switch_param (must match existing switch_param), time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same ChemotaxA experiment
    def chemotax_heatmap(self, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        chem_num = len(self.exp_type[0].chemotax_params)
        move_axis = np.zeros((chem_num + 1, move_num + 1))
        chem_axis = np.zeros((move_num + 1, chem_num + 1))
        for ch_index in range(chem_num+1):
            move_axis[ch_index] = np.asarray(self.exp_type[0].move_params+[31.6])
        for m_index in range(move_num+1):
            chem_axis[m_index] = np.asarray([chem-0.05 for chem in self.exp_type[0].chemotax_params]+[0.95])
        chem_axis = np.transpose(chem_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, chem_num))
        for move_index in range(move_num):
            for chem_index in range(chem_num):
                move = self.exp_type[0].move_params[move_index]
                chem = self.exp_type[0].chemotax_params[chem_index]
                if self.adap_rate((move, chem)) is None:
                    adap_rate = 0
                else:
                    adap_rate = self.adap_rate((move, chem))[0][time_type][LD-1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[move_index][chem_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        adap_rates = np.transpose(adap_rates)
        im = ax.pcolormesh(move_axis, chem_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.axvline(self.exp_type[0].DeathRate[0], color='black', linestyle="--")
        ax.text(0.5 * self.exp_type[0].DeathRate[0], 0.07, "$\\delta$", fontsize=12, color='black')
        ax.axhline(0.5, color='black', linestyle="--")
        ax.text(1.4*self.exp_type[0].move_params[0], 0.52, "positive  chemotaxis $\\uparrow$", fontsize=12, color='black')
        ax.text(1.4 * self.exp_type[0].move_params[0], 0.45, "negative chemotaxis $\\downarrow$", fontsize=12, color='black')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            if labels:
                cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
            ax.set_ylabel("chemotaxis p", fontsize=12)
        return [min, max, im]

    # make the heatmap plot of adaptation rat of given time_type, LD
    # output: [min_adap_rate, max_adap_rate, im]
    # input: switch_param (must match existing switch_param), time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same ResCostA or ResCostDeadlyA experiment
    def res_cost_heatmap(self, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        cost_num = len(self.exp_type[0].cost_params)
        move_axis = np.zeros((cost_num + 1, move_num + 1))
        cost_axis = np.zeros((move_num + 1, cost_num + 1))
        for c_index in range(cost_num+1):
            move_axis[c_index] = np.asarray(self.exp_type[0].move_params+[31.6])
        for m_index in range(move_num+1):
            cost_axis[m_index] = np.asarray(self.exp_type[0].cost_params+[1])
        cost_axis = np.transpose(cost_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, cost_num))
        for move_index in range(move_num):
            for cost_index in range(cost_num):
                move = self.exp_type[0].move_params[move_index]
                cost = self.exp_type[0].cost_params[cost_index]
                if self.adap_rate((move, cost)) is None:
                    adap_rate = 0
                else:
                    adap_rate = self.adap_rate((move, cost))[0][time_type][LD-1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[move_index][cost_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        adap_rates = np.transpose(adap_rates)
        im = ax.pcolormesh(move_axis, cost_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(self.exp_type[0].DeathRate[0], color='black', linestyle="--")
        ax.text(0.5 * self.exp_type[0].DeathRate[0], 1e-4, "$\\delta$", fontsize=12, color='black')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            cbar.set_label("adaptation rate $a$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
            ax.set_ylabel("resistance cost $c$", fontsize=12)
        return [min, max, im]

    # make the heatmap plot of adaptation rat of given time_type, LD
    # output: [min_adap_rate, max_adap_rate, im]
    # input: switch_param (must match existing switch_param), time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: firstly input of the same HGTA experiment or the same HGTDeadlyA experiment; then can include standard
    def hgt_heatmap(self, LD, time_type, fig, ax, min_max, colorbar, labels, line):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        hgt_num = len(self.exp_type[0].hgt_params)
        if isinstance(self.exp_type[-1], StandardA):
            hgt_num += 1
        move_axis = np.zeros((hgt_num + 1, move_num + 1))
        hgt_axis = np.zeros((move_num + 1, hgt_num + 1))
        for h_index in range(hgt_num+1):
            move_axis[h_index] = np.asarray(self.exp_type[0].move_params+[31.6])
        for m_index in range(move_num+1):
            if isinstance(self.exp_type[-1], StandardA):
                hgt_axis[m_index] = np.asarray([3.16e-8] + [hgt for hgt in self.exp_type[0].hgt_params] + [3.16e-3])
            else:
                hgt_axis[m_index] = np.asarray([hgt for hgt in self.exp_type[0].hgt_params]+[3.16e-3])
        hgt_axis = np.transpose(hgt_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, hgt_num))
        for move_index in range(move_num):
            for hgt_index in range(hgt_num):
                params = []
                if isinstance(self.exp_type[-1], StandardA):
                    if hgt_index == 0:
                        move = self.exp_type[0].move_params[move_index]
                        death = self.exp_type[0].DeathRate[0]
                        for exp in self.exp_type:
                            if isinstance(exp, StandardA):
                                params.append((move, death))
                            else:
                                params.append(None)
                    else:
                        move = self.exp_type[0].move_params[move_index]
                        hgt = self.exp_type[0].hgt_params[hgt_index-1]
                        for exp in self.exp_type:
                            if isinstance(exp, StandardA):
                                params.append(None)
                            else:
                                params.append((move, hgt))
                else:
                    move = self.exp_type[0].move_params[move_index]
                    hgt = self.exp_type[0].hgt_params[hgt_index]
                    for exp in self.exp_type:
                            params.append((move, hgt))
                if self.adap_rate(params) is None:
                    adap_rate = 0
                else:
                    adap_rate = self.adap_rate(params)[0][time_type][LD-1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[move_index][hgt_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        adap_rates = np.transpose(adap_rates)
        im = ax.pcolormesh(move_axis, hgt_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(self.exp_type[0].DeathRate[0], color='black', linestyle="--")
        ax.text(0.5 * self.exp_type[0].DeathRate[0], 2e-7, "$\\delta$", fontsize=12, color='black')
        if self.exp_type[0].DeathRate[0]==3e-1:
            ax.axvline(3.16, color='black', linestyle="--")
            ax.text(0.25 * 3.16, 2e-7, "$\\nu_c$", fontsize=12, color='black')
        if isinstance(self.exp_type[-1], StandardA):
            ax.axhline(1e-7, color='black', linestyle="-")
            ax.text(2e-5, 4e-8, "NO HGT", fontsize=12, color='white')
            ax.minorticks_off()
            yvalues = [1e-7, 1e-6, 1e-5, 1e-4]
            xvalues = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10]
            yticks = []
            xticks = []
            for value in yvalues:
                yticks = yticks+[value*(i+2) for i in range(8)]
            yticks = yticks + [2e-3, 3e-3]
            for value in xvalues:
                xticks = xticks + [value * (i + 2) for i in range(8)]
            yticks = yticks + [20, 30]
            ax.set_yticks(yticks, minor=True)
            ax.set_xticks(xticks, minor=True)
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if line:
            mov1=1e-5
            mov2=31.6
            hgt1=1e-4
            hgt2=3.16e-7
            b=np.log(hgt2/hgt1)/(mov2-mov1)
            a=hgt1*np.exp(-mov1*b)
            x_axis = np.linspace(mov1,mov2,num=100)
            y_axis = np.asarray([a*np.exp(b*mov) for mov in x_axis])
            ax.plot(x_axis, y_axis, color='grey', linewidth=3)
            ax.text(3.16e-4, 2e-4, "PHYSICAL\nSYSTEM", color="grey", fontsize=12)
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
            ax.set_ylabel("HGT rate $\\mu_{HGT}$ (1/h)", fontsize=12)
        return [min, max, im]

    # make the heatmap plot of adaptation rate of given time_type, LD
    # output: [min_adap_rate, max_adap_rate, im]
    # input: switch_param (must match existing switch_param), time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same StressDeathB experiment
    def stressdeath_heatmap(self, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        sdeath_num = len(self.exp_type[0].sdeath_params)
        move_axis = np.zeros((sdeath_num + 1, move_num + 1))
        sdeath_axis = np.zeros((move_num + 1, sdeath_num + 1))
        for sd_index in range(sdeath_num+1):
            move_axis[sd_index] = np.asarray(self.exp_type[0].move_params+[31.6])
        for m_index in range(move_num+1):
            sdeath_axis[m_index] = np.asarray([(sdeath+0.1)/0.1-0.5 for sdeath in self.exp_type[0].sdeath_params]+[9.5])
        sdeath_axis = np.transpose(sdeath_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, sdeath_num))
        for move_index in range(move_num):
            for sdeath_index in range(sdeath_num):
                move = self.exp_type[0].move_params[move_index]
                sdeath = self.exp_type[0].sdeath_params[sdeath_index]
                if self.adap_rate((move, sdeath)) is None:
                    adap_rate = 0
                else:
                    adap_rate = self.adap_rate((move, sdeath))[0][time_type][LD-1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[move_index][sdeath_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        adap_rates = np.transpose(adap_rates)
        im = ax.pcolormesh(move_axis, sdeath_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.axvline(self.exp_type[0].DeathRate[0], color='black', linestyle="--")
        ax.text(0.5 * self.exp_type[0].DeathRate[0], 0.7, "$\\delta$", fontsize=12, color='black')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            if labels:
                cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)", fontsize=12)
            ax.set_ylabel("drug-induced death ratio $\\sigma$", fontsize=12)
        return [min, max, im]

    # make the heatmap plot of adaptation rate of given time_type, LD, and switch parameter
    # output: [min_adap_rate, max_adap_rate, im]
    # input: switch_param (must match existing switch_param), time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same MotSwitchA experiment
    def mot_switch_heatmap(self, switch_param, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        move0_axis = np.zeros((move_num + 1, move_num + 1))
        for m_index in range(move_num+1):
            move0_axis[m_index] = np.asarray(self.exp_type[0].move_params+[31.6])
        move1_axis = copy.deepcopy(move0_axis)
        move1_axis = np.transpose(move1_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, move_num))
        for move0_index in range(move_num):
            for move1_index in range(move_num):
                move0 = self.exp_type[0].move_params[move0_index]
                move1 = self.exp_type[0].move_params[move1_index]
                if move0 < move1:
                    if self.adap_rate((move1, move0, switch_param)) is None:
                        adap_rate = 0
                    else:
                        adap_rate = self.adap_rate((move1, move0, switch_param))[0][time_type][LD-1]
                        if np.isnan(adap_rate):
                            adap_rate = 0
                    adap_rates[move0_index][move1_index] = adap_rate
                else:
                    if self.adap_rate((move0, move1, switch_param)) is None:
                        adap_rate = 0
                    else:
                        adap_rate = self.adap_rate((move0, move1, switch_param))[0][time_type][LD - 1]
                        if np.isnan(adap_rate):
                            adap_rate = 0
                    adap_rates[move0_index][move1_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min / 2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        im = ax.pcolormesh(move0_axis, move1_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.set_yscale("log")
        # add guiding lines nu=delta, nu1+nu2=delta
        ax.vlines(self.exp_type[0].DeathRate[0], self.exp_type[0].DeathRate[0], 31.6, color='black', linestyle="--")
        #ax.text(0.3*self.exp_type[0].DeathRate[0], 1.4*self.exp_type[0].move_params[0], "$\\delta$", fontsize=12, color='black')
        ax.hlines(self.exp_type[0].DeathRate[0], self.exp_type[0].DeathRate[0], 31.6, color='black', linestyle="--")
        #ax.text(1.4*self.exp_type[0].move_params[0], 0.3 * self.exp_type[0].DeathRate[0], "$\\delta$", fontsize=12, color='black')
        nu1_vals = np.logspace(-5, np.log10(2*self.exp_type[0].DeathRate[0]), num=100)
        nu2_vals = [2*self.exp_type[0].DeathRate[0]-nu1 for nu1 in nu1_vals]
        ax.plot(nu1_vals, nu2_vals, color="black", linestyle='--')
        # add guiding lines at low switching
        if switch_param < self.exp_type[0].DeathRate[0]:
            for exponent in [-5,-4,-3,-2,-1,0,1]:
                ax.vlines(10**exponent, 10**exponent, 31.6, color='white', linestyle="--")
                ax.hlines(10**exponent, 10**exponent, 31.6, color='white', linestyle="--")
        else:
            for exponent in [-5, -4, -3, -2, -1,0,1]:
                nu1_vals = np.logspace(-5, np.log10(2*(10**exponent)), num=100)
                nu2_vals = [2*(10**exponent)-nu1 for nu1 in nu1_vals]
                ax.plot(nu1_vals, nu2_vals, color="white", linestyle='--')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            if labels:
                cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("motility $\\nu_1$ (1/h)", fontsize=12)
            ax.set_ylabel("motility $\\nu_2$ (1/h)", fontsize=12)
            ax.set_title("switching rate $s=$"+str(switch_param))
        return [min, max, im]

    # make the heatmap plot of adaptation rate of given time_type, LD, and switch_dens parameter, death parameter (must match existing param)
    # output: [min_adap_rate, max_adap_rate, im]
    # input: dens_param, switch_dens, time_type (0=T, 1=D), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same DensSwitchA experiment
    def dens_switch_heatmap(self, death, dens_switch, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].move_params)
        move0_axis = np.zeros((move_num + 1, move_num + 1))
        for m_index in range(move_num + 1):
            move0_axis[m_index] = np.asarray(self.exp_type[0].move_params + [31.6])
        move1_axis = copy.deepcopy(move0_axis)
        move1_axis = np.transpose(move1_axis)
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((move_num, move_num))
        for move0_index in range(move_num):
            for move1_index in range(move_num):
                move0 = self.exp_type[0].move_params[move0_index]
                move1 = self.exp_type[0].move_params[move1_index]
                if self.adap_rate((death, dens_switch, move0, move1)) is None:
                    adap_rate = 0
                else:
                    adap_rate = self.adap_rate((death, dens_switch, move0, move1))[0][time_type][LD - 1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[move1_index][move0_index] = adap_rate
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min / 2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min / 2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        im = ax.pcolormesh(move0_axis, move1_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.set_yscale("log")
        # add guiding lines nu=delta, nu1+nu2=delta
        ax.vlines(death, 1e-5, 31.6, color='black', linestyle="--")
        ax.text(0.3 * death, 1.4 * 1e-5, "$\\delta$", fontsize=12, color='black')
        ax.hlines(death, 1e-5, 31.6, color='black', linestyle="--")
        ax.text(1.4 * 1e-5, 0.3 * death, "$\\delta$", fontsize=12, color='black')
        ax.plot([1e-5, 31.6], [1e-5, 31.6], linestyle=":", color='black')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            if labels:
                cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("motility $\\nu_1$ (1/h)", fontsize=12)
            ax.set_ylabel("motility $\\nu_2$ (1/h)", fontsize=12)
            ax.set_title("switching threshold $S=$" + str(dens_switch))
        return [min, max, im]

    # swarm adap plot
    # output: []
    # input: part(standard=0,swarm=1), death, swarm_density, LD, min_max, ax, labels, background
    # note: assume all input data EITHER StandardA experiment OR all input data SwarmDeadlyA experiment
    def swarm_adap_plot(self, part, death, swarm_density, LD, min_max, ax, labels, background):
        swarm_deadly = SwarmDeadlyA
        move_num = len(swarm_deadly.swarm_move_params)
        adap_rates = np.zeros((2, move_num))
        crit_mot = None

        # Standard A:
        if part == 0:
            # Swimming adap rates
            for move_index in range(move_num):
                move = swarm_deadly.swarm_move_params[move_index]
                # get adaptation rate
                adap_rate = self.adap_rate((move, death))
                if adap_rate is None:
                    adap_rate = 'NaN'
                else:
                    adap_rate = adap_rate[0][0][LD - 1]
                    if np.isnan(adap_rate):
                        adap_rate = 'NaN'
                adap_rates[0][move_index] = adap_rate
                # get critical motility
                if crit_mot is None:
                    if adap_rate == 'NaN':
                        crit_mot = move
                else:
                    if adap_rate != 'NaN':
                        crit_mot = None
            # Stationary adap rates
            adap_rate = self.adap_rate((1e-3, death))[0][0][LD-1]
            adap_rates[1] = np.asarray([adap_rate for i in range(move_num)])
            ax.plot(swarm_deadly.swarm_move_params, adap_rates[0], linestyle='--', marker='o', label="faster only", color="#D4AC0D")
            ax.plot(swarm_deadly.swarm_move_params, adap_rates[1], linestyle='--', marker='o', label="slower only", color="#2E86C1")
        # SwarmDeadlyA
        else:
            for move_index in range(move_num):
                move = swarm_deadly.swarm_move_params[move_index]
                # get adaptation rate
                adap_rate = self.adap_rate((move, swarm_density, death))
                if adap_rate is None:
                    adap_rateT = 'NaN'
                    adap_rateP = 'NaN'
                else:
                    adap_rateT = adap_rate[0][0][LD-1]
                    adap_rateP = adap_rate[0][2][1]
                    if np.isnan(adap_rateT):
                        adap_rateT = 'NaN'
                    if np.isnan(adap_rateP):
                        adap_rateP = 'NaN'
                adap_rates[0][move_index] = adap_rateT
                adap_rates[1][move_index] = adap_rateP
            ax.plot(swarm_deadly.swarm_move_params, adap_rates[1], marker='o', label="swarming tolerance", color="tab:green")
            ax.plot(swarm_deadly.swarm_move_params, adap_rates[0], marker='o', label="genetic resistance", color="tab:red")

        # set common features of the plot
        if min_max is None:
            max = adap_rates.max()
            min = np.nanmin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylim([min,max])
        ax.set_xlim(swarm_deadly.swarm_move_params[0]*0.5, swarm_deadly.swarm_move_params[-1]*1.5)
        # make the background
        if background:
            if crit_mot is None:
                ax.axvspan(swarm_deadly.swarm_move_params[0] / 5, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, swarm_deadly.swarm_move_params[-1] * 5, facecolor="tab:orange", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.6 * death, 6e-4, "$\\delta$", fontsize=12, color='black')
                #ax.text(swarm_deadly.swarm_move_params[0] / 2.5, 6e-3, "Low Motility", fontsize=15, color='green', rotation=90)
                #ax.text(swarm_deadly.swarm_move_params[-1]*2.5, 6e-3, "High Motility", fontsize=15, color='orange', rotation=90)
                ax.text(swarm_deadly.swarm_move_params[0], 3e-2, "LM", fontsize=15, color='green')
                ax.text(swarm_deadly.swarm_move_params[-1]/10, 3e-2, "HM", fontsize=15, color='orange')
            else:
                ax.axvspan(swarm_deadly.swarm_move_params[0] / 5, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, crit_mot, facecolor="tab:orange", alpha=0.3)
                ax.axvspan(crit_mot, swarm_deadly.swarm_move_params[-1] * 5, facecolor="tab:red", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.6 * death, 6e-4, "$\\delta$", fontsize=12, color='black')
                ax.axvline(crit_mot, color='black', linestyle="--")
                ax.text(1.1 * crit_mot, 6e-4, "$\\nu_c$", fontsize=12)
                #ax.text(swarm_deadly.swarm_move_params[0] / 2.5, 6e-3, "Low Motility", fontsize=15, color='green', rotation=90)
                #ax.text(2 * death, 6e-3, "High Motility", fontsize=15, color='orange', rotation=90)
                #ax.text(swarm_deadly.swarm_move_params[-1]*2.5, 6e-3, "Deadly Motility", fontsize=15, color='red', rotation=90)
                ax.text(swarm_deadly.swarm_move_params[0], 3e-2, "LM", fontsize=15, color='green')
                ax.text(1.5 * death, 3e-2, "HM", fontsize=15, color='orange')
                ax.text(swarm_deadly.swarm_move_params[-1]/2.5, 3e-2, "DM", fontsize=15, color='red')
        # plot the labels
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)")
            ax.set_ylabel("adaptation rates $a_2$, $a_S$ (1/h)")

    # make the heatmap plot of adaptation rate of given time_type, LD, and swarm parameter
    # output: [min_adap_rate, max_adap_rate, im]
    # input: death rate, time_type (0=T, 1=D, 2=PA), LD, fig, ax, min_max, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same SwarmDeadlyA experiment
    def swarm_adap_heatmap(self, death, LD, time_type, fig, ax, min_max, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].swarm_move_params)
        density_num = len(self.exp_type[0].swarm_density_params)
        move_axis = np.asarray(self.exp_type[0].swarm_move_params+[31.6])
        density_axis = np.asarray(self.exp_type[0].swarm_density_params+[1e5])
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        adap_rates = np.zeros((density_num, move_num))
        for move_index in range(move_num):
            for dense_index in range(density_num):
                move = self.exp_type[0].swarm_move_params[move_index]
                dense = self.exp_type[0].swarm_density_params[dense_index]
                adap_rate = self.adap_rate((move, dense, death))
                if adap_rate is None:
                    adap_rate = 0
                else:
                    if time_type == 2:
                        adap_rate = adap_rate[0][time_type][LD]
                    else:
                        adap_rate = adap_rate[0][time_type][LD-1]
                    if np.isnan(adap_rate):
                        adap_rate = 0
                adap_rates[dense_index][move_index] = adap_rate
        print(adap_rates)
        if min_max is None:
            max = adap_rates.max()
            if np.amin(adap_rates) == 0:
                min = np.amin(np.array(adap_rates)[adap_rates != 0])
                adap_rates = np.where(adap_rates == 0, min / 2, adap_rates)
            else:
                min = np.amin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
            if np.amin(adap_rates) == 0:
                adap_rates = np.where(adap_rates == 0, min/2, adap_rates)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        im = ax.pcolormesh(move_axis, density_axis, adap_rates, cmap=my_cmap, norm=LogNorm(vmin=min, vmax=max))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(death, color='black', linestyle="--")
        ax.text(0.3*death, 1.4*self.exp_type[0].swarm_move_params[0], "$\\delta$", fontsize=12, color='black')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            if labels:
                if time_type == 2:
                    cbar.set_label("swarming rate $a_S$ (1/h)", fontsize=12, rotation=90)
                else:
                    cbar.set_label("adaptation rate $a_2$ (1/h)", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("swarming motility $\\nu_S$ (1/h)", fontsize=12)
            ax.set_ylabel("swarming density $K_S$", fontsize=12)
        return [min, max, im]

    # make the heatmap plot of expected phen. tolerance occured before gen. res. of given time_type, LD, and swarm parameter
    # output: [min_adap_rate, max_adap_rate, im]
    # input: death rate, LD, swarming population, fig, ax, colorbar (True/False), labels(True,False)
    # note: assume all input data are of the same SwarmDeadlyA experiment
    def swarm_tol_heatmap(self, death, LD, pop, fig, ax, colorbar, labels):
        # prepare x, y axes
        move_num = len(self.exp_type[0].swarm_move_params)
        density_num = len(self.exp_type[0].swarm_density_params)
        move_axis = np.asarray(self.exp_type[0].swarm_move_params+[31.6])
        density_axis = np.asarray(self.exp_type[0].swarm_density_params+[1e5])
        # prepare adaptation rates, and find minimal non-zero and maximal adaptation rate
        swarm_tol = np.zeros((density_num, move_num))
        for move_index in range(move_num):
            for dense_index in range(density_num):
                move = self.exp_type[0].swarm_move_params[move_index]
                dense = self.exp_type[0].swarm_density_params[dense_index]
                swarm_tol[dense_index][move_index] = self.swarm_tolerance((move, dense, death), pop)[0][LD-1]
        print(swarm_tol)
        # prepare colormap which has a red color for no adaptation rate
        my_cmap = cm.get_cmap('viridis').copy()
        my_cmap.set_under(color='red')
        # make plot and add features
        im = ax.pcolormesh(move_axis, density_axis, swarm_tol, cmap=my_cmap, vmin=0, vmax=self.exp_type[0].GridBound-1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axvline(death, color='black', linestyle="--")
        ax.text(1.5*death, 1.2*self.exp_type[0].swarm_density_params[0], "$\\delta$", fontsize=12, color='black')
        if death==3e-1:
            ax.axvline(3, color='black', linestyle="--")
            ax.text(1.2, 1.2*self.exp_type[0].swarm_density_params[0], "$\\nu_c$", fontsize=12, color='black')
        # plot threshold lines
        motility1 = np.logspace(np.log10(self.exp_type[0].swarm_move_params[0]), np.log10(death), num=100)
        KS1 = motility1/death*self.exp_type[0].CarryingCapacity*(1-death/self.exp_type[0].BirthRate[0])-motility1**2*self.exp_type[0].CarryingCapacity/(death*self.exp_type[0].BirthRate[0])
        motility2 = np.logspace(np.log10(death), np.log10(31.6), num=100)
        KS2 = [self.exp_type[0].CarryingCapacity*(1-(death+self.exp_type[0].Move[0])/self.exp_type[0].BirthRate[0]) for i in range(100)]
        ax.plot(motility1, KS1, color='red', linestyle='--')
        ax.plot(motility2, KS2, color='red', linestyle='--')
        # set colorbar and labels
        if colorbar:
            cbar = fig.colorbar(im)
            cbar.set_label("swarming tolerance", fontsize=12, rotation=90)
        if labels:
            ax.set_xlabel("swarming motility $\\nu_S$ (1/h)", fontsize=12)
            ax.set_ylabel("swarming density $K_S$", fontsize=12)
        return im

    # make the plot of swarming radius vs time for a given innoculation density and swarming motility
    # output: im, time_axis, radius_axis
    # input: init_cell, swarm_move, averaged (True,False), fig, ax, labels(True,False)
    # notice: init_cell number = [1e2, 2.15e2, 4.64e2, 1e3, 2.15e3, 4.64e3, 1e4, 2.15e4, 4.64e4]
    # note: assume all input data are of the same SwarmLagA experiment
    def swarm_radius_plot(self, init_cell, swarm_move, averaged, ax, labels):
        # prepare values for plotting axes
        time_axis = [0]
        radius_axis = [0]
        time = 0
        if averaged[0]:
            wait_times = self.waiting_time((swarm_move, init_cell))
            print(wait_times[0])
            for LD in range(len(wait_times[0][2])):
                wait_time = wait_times[0][2][LD]
                if not np.isnan(wait_time):
                    time += wait_time
                    time_axis.append(time)
                    radius_axis.append(LD)
        else:
            task_id = self.exp_type[0].IDtoTaskID((swarm_move, init_cell), averaged[2])
            times = self.read_times(averaged[1],task_id)
            print(times)
            if times is not None:
                times = times[2][0]
                for LD in range(len(times)):
                    time_axis.append(times[LD])
                    radius_axis.append(LD)
        im = ax.plot(time_axis, radius_axis, color="black", linewidth=2)
        # set colorbar and labels
        if labels:
            ax.set_xlabel("time (h)", fontsize=12)
            ax.set_ylabel("swarming radius", fontsize=12)
        return [im, time_axis, radius_axis]

    # make motility-adaptation rate plot for special
    # output: [min_adap_rate, max_adap_rate]
    # input: special_type (0=res. cost, 1=sbirth, 2=sdeath), time_type (0=T, 1=D), LD,
    #        analytics(t/f, poin_num, precision), ax, min_max, labels (t/f), legend_name (string), background (t/f)
    # note: assume all input data are of the SpecialA or SpecialB experiment
    def special_motility_adap_plot(self, special_type, LD, time_type, analytics, ax, min_max, labels, legend_name, background):
        # prepare adaptation rates, and find minimal and maximal adaptation rate
        move_num = len(self.exp_type[0].move_params)
        adap_rates = np.zeros(move_num)
        crit_mot = None
        for move_index in range(move_num):
            move = self.exp_type[0].move_params[move_index]
            # get adaptation rate
            adap_rate = self.adap_rate((move, special_type))
            if adap_rate is None:
                adap_rate = 'NaN'
            else:
                adap_rate = adap_rate[0][time_type][LD - 1]
                if np.isnan(adap_rate):
                    adap_rate = 'NaN'
            adap_rates[move_index] = adap_rate
            # get critical motility
            if crit_mot is None:
                if adap_rate == 'NaN':
                    crit_mot = move
            else:
                if adap_rate != 'NaN':
                    crit_mot = None
        if min_max is None:
            max = adap_rates.max()
            min = np.nanmin(adap_rates)
        else:
            min = min_max[0]
            max = min_max[1]
        # make the plot
        ax.scatter(self.exp_type[0].move_params, adap_rates, label=legend_name)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(self.exp_type[0].move_params[0] / 10, self.exp_type[0].move_params[-1] * 10)
        # make the background
        if background:
            death = self.exp_type[0].DeathRate[0]
            if crit_mot is None:
                ax.axvspan(self.exp_type[0].move_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, self.exp_type[0].move_params[-1] * 10, facecolor="tab:orange", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.1*min, "$\\delta$", fontsize=12, color='black')
                ax.text(self.exp_type[0].move_params[0] / 5, 1.5 * min, "Low Motility", fontsize=12, color='green', rotation=90)
                ax.text(2 * death, 1.5 * min, "High Motility", fontsize=12, color='orange', rotation=90)
            else:
                ax.axvspan(self.exp_type[0].move_params[0] / 10, death, facecolor="tab:green", alpha=0.3)
                ax.axvspan(death, crit_mot, facecolor="tab:orange", alpha=0.3)
                ax.axvspan(crit_mot, self.exp_type[0].move_params[-1] * 10, facecolor="tab:red", alpha=0.3)
                ax.axvline(death, color='black', linestyle="--")
                ax.text(0.4 * death, 1.1*min, "$\\delta$", fontsize=12, color='black')
                ax.axvline(crit_mot, color='black', linestyle="-.")
                ax.text(1.2 * crit_mot, 1.1*min, "$\\nu_c$", fontsize=12)
                ax.text(self.exp_type[0].move_params[0] / 5, 1.5 * min, "Low Motility", fontsize=12, color='green', rotation=90)
                ax.text(2 * death, 1.5 * min, "High Motility", fontsize=12, color='orange', rotation=90)
                ax.text(5 * crit_mot, 1.5 * min, "Deadly Motility", fontsize=12, color='red', rotation=90)
        # plot the labels
        if labels:
            ax.set_xlabel("motility $\\nu$ (1/h)")
            ax.set_ylabel("adaptation $a_2$ (1/h)")
        # make and plot the analytics
        if analytics[0]:
            point_num = analytics[1]
            move_params = np.logspace(np.log10(self.exp_type[0].move_params[0]/2), np.log10(2*self.exp_type[0].move_params[-1]), num=point_num)
            adap_params = np.zeros(point_num)
            l_prop = sc.LatProp()
            l_prop.LD = LD
            p_prop = [sc.PopProp()]
            if special_type == 0:
                p_prop[0].ResistCost = self.exp_type[0].special_params[0]
            elif special_type == 1:
                p_prop[0].StressBirthRate = self.exp_type[0].special_params[1]
            elif special_type == 2:
                p_prop[0].StressDeathRate = self.exp_type[0].special_params[2]
                if isinstance(self.exp_type[0], SpecialB):
                    p_prop[0].StressBirthRate = self.exp_type[0].special_params[3]
            for i in range(point_num):
                p_prop[0].Move = move_params[i]
                stair = sc.Staircase(l_prop, p_prop)
                adap_params[i] = stair.retarded_adap_rate(analytics[2], 1e4)
            ax.plot(move_params, adap_params, label="")
        return[min, max]

    # make RC/SB/SD-adaptation rate plot
    # output: [min_adap_rate, max_adap_rate]
    # input: time_type (0=T, 1=D), LD, analytics (t/f, poin_num, precision), ax, labels (t/f)
    # note: accepts the following combinations of experiments
    # 1) ResistCostA,B (mixed), 2) StressBirthA (only), 3) StressDeathA (only), 4) StressBirthDeathA
    def special_adap_plot(self, LD, time_type, analytics, ax, labels):
        # get general staircase information (that does not change)
        vary_params = None
        all_rc_A = True
        rc_B_index = None
        type = 0  # 0=RC, 1=SB, 2=SD, 3=SBD
        if isinstance(self.exp_type[0], ResistCostA) or isinstance(self.exp_type[0], ResistCostB):
            for i in range(len(self.exp_type)):
                if isinstance(self.exp_type[i], ResistCostB):
                    all_rc_A = False
                    rc_B_index = i
            if all_rc_A:
                vary_params = self.exp_type[0].cost_params
            else:
                vary_params = self.exp_type[rc_B_index].cost_params
            type = 0
        elif isinstance(self.exp_type[0], StressBirthA):
            vary_params = self.exp_type[0].sbirth_params
            type = 1
        elif isinstance(self.exp_type[0], StressDeathA):
            vary_params = self.exp_type[0].sdeath_params
            type = 2
        elif isinstance(self.exp_type[0], StressBirthDeathA):
            vary_params = self.exp_type[0].sdeath_params
            type = 3
        else:
            "Warning! Wrong experiment choices."
        # prepare adaptation rates, and find minimal and maximal adaptation rate
        params = None
        vary_num = len(vary_params)
        mov_num = len(self.exp_type[0].move_params)
        adap_rates = np.zeros((mov_num, vary_num))
        crit_mot = [None for m_i in range(mov_num)]
        for move_index in range(mov_num):
            for vary_index in range(vary_num):
                move = self.exp_type[0].move_params[move_index]
                vary = vary_params[vary_index]
                # get correct params if all_rc_A is false
                if (not all_rc_A) and (vary == 6.3e-1 or vary == 7.94e-1):
                    params = []
                    for i in range(len(self.exp_type)):
                        if isinstance(self.exp_type[i], ResistCostA):
                            params.append(None)
                        else:
                            params.append((move, vary))
                else:
                    params = (move, vary)
                # get adaptation rate
                adap_rate = self.adap_rate(params)
                if adap_rate is None:
                    adap_rate = 'NaN'
                else:
                    adap_rate = adap_rate[0][time_type][LD - 1]
                    if np.isnan(adap_rate):
                        adap_rate = 'NaN'
                adap_rates[move_index][vary_index] = adap_rate
                # get critical motility
                if crit_mot[move_index] is None:
                    if adap_rate == 'NaN':
                        crit_mot[move_index] = move
                else:
                    if adap_rate != 'NaN':
                        crit_mot[move_index] = None
        # make the plot
        if type == 3:
            vary_params = [(vary+0.1)*10 for vary in vary_params]
        for move_index in range(mov_num):
            if self.exp_type[0].move_params[move_index]<self.exp_type[0].DeathRate[0]:
                ax.scatter(vary_params, adap_rates[move_index], label="$\\nu \\ll \\delta$")
            else:
                ax.scatter(vary_params, adap_rates[move_index], label="$\\nu \\gg \\delta$")
        ax.set_yscale("log")
        ax.set_xscale("log")
        # plot the labels
        if labels:
            if type == 0:
                ax.set_xlabel("resistance cost $c$")
            elif type == 1:
                ax.set_xlabel("birth under stress $r_S$ (1/h)")
            elif type == 2 or type == 3:
                ax.set_xlabel("stress-induced death $\\delta_S$ (1/h)")
            ax.set_ylabel("adaptation $a_2$ (1/h)")
        ax.legend(loc='lower left')
        # analytics
        if analytics[0]:
            # prepare x-axis
            point_num = analytics[1]
            varying_params = None
            if type == 0:
                print("RESISTANCE COST")
                if all_rc_A:
                    varying_params = np.logspace(-4, -0.06, num=point_num)
                else:
                    varying_params = np.logspace(-4, -0.06, num=point_num)
            elif type == 1:
                print("STRESS BIRTH")
                varying_params = np.logspace(-3, -0.06, num=point_num)
            elif type == 2:
                print("STRESS DEATH")
                varying_params = np.logspace(-1, 0, num=point_num)
            elif type == 3:
                print("STRESS BIRTH DEATH")
                varying_params = np.logspace(1, 2, num=point_num)
            # prepare both y-axes values
            l_prop = sc.LatProp()
            l_prop.LD = LD
            p_prop = [sc.PopProp()]
            adap_params = np.zeros((mov_num, point_num))
            for i in range(point_num):
                if type == 0:
                    p_prop[0].ResistCost = varying_params[i]
                elif type == 1:
                    p_prop[0].StressBirthRate = varying_params[i]
                elif type == 2:
                    p_prop[0].StressDeathRate = varying_params[i]
                elif type == 3:
                    p_prop[0].StressBirthRate = 1
                    p_prop[0].StressDeathRate = varying_params[i]/10-0.1
                for j in range(mov_num):
                    p_prop[0].Move = self.exp_type[0].move_params[j]
                    stair = sc.Staircase(l_prop, p_prop)
                    adap_params[j][i] = stair.retarded_adap_rate(analytics[2], 1e4)
            ax.plot(varying_params, adap_params[0], label="")
            ax.plot(varying_params, adap_params[1], label="")
        return [min, max]

    # plots the staircase plot into given axis
    # input: params, run_index(=0,...,repeat), ax(axes to plot to), labels(=true,false)
    #        sim_id, evol_pop (just evolved), plot_pop (to be plotted), LD, time_type(='T','D')
    def staircase_plot(self, params, run_index, sim_id, evol_pop, plot_pop, LD, time_type, ax, labels):
        # get important data
        task_id = self.exp_type[sim_id].IDtoTaskID(params, run_index)
        [Time, GenSpaceTot] = self.read_gen_space_tot(sim_id, task_id, evol_pop, LD, time_type)
        # make the plot
        im = ax.imshow(GenSpaceTot[plot_pop], vmin=0, vmax=self.exp_type[sim_id].CarryingCapacity, origin='lower', cmap='Greys')
        for j in range(self.exp_type[sim_id].GridBound-1):
            ax.hlines(j + 1 - 0.5, j + 1 - 0.5, j + 2 - 0.5, colors='black')
            ax.vlines(j + 1 - 0.5, j - 0.5, j + 1 - 0.5, colors='black')
        # make lables
        if not labels:
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.text(0, self.exp_type[sim_id].GenBound, '$L_D=$' + str(LD) + ", t=" + str(Time))
            ax.set_xlabel('space x')
            ax.set_ylabel('genotype g')

    # calculate space_tot from gen_space_tot
    def space_tot(self, gen_space_tot, sim_id):
        space_tot = np.zeros((self.exp_type[sim_id].PopulationNumber, self.exp_type[sim_id].GridBound))
        for pop in range(self.exp_type[sim_id].PopulationNumber):
            for gen in range(self.exp_type[sim_id].GenBound):
                for pos in range(self.exp_type[sim_id].GridBound):
                    space_tot[pop][pos] += gen_space_tot[pop][gen][pos]
        return space_tot

    # finds the swarming number S for a given swarm_pop and GenSpaceTot
    def S(self, gen_space_tot, swarm_pop, swarm_density, sim_id):
        if self.exp_type[sim_id].ConsiderSwarm[swarm_pop]:
            S = 0
            space_tot = self.space_tot(gen_space_tot, sim_id)[swarm_pop]
            while S < self.exp_type[sim_id].GridBound and space_tot[S] >= swarm_density * 0.99:
                S += 1
            return S
        else:
            return None

    # plots the stable wild-type plot into given axis, colors indicate relative abundance
    # input: params, run_index(=0,...,repeat), ax(axes to plot to), labels(=true,false)
    #        sim_id, evol_pop (just evolved), show_r_net, LD, time_type(='T','D'), ax, labels, analytics
    def wild_type_plot(self, params, run_index, sim_id, evol_pop, show_r_net, LD, time_type, ax, x_labels, pop_labels, r_labels, analytics):
        # get important data
        task_id = self.exp_type[sim_id].IDtoTaskID(params, run_index)
        if self.read_gen_space_tot(sim_id, task_id, evol_pop, LD, time_type) is None:
            [Time,GenSpaceTot] = [np.inf, np.zeros((self.exp_type[sim_id].PopulationNumber, self.exp_type[sim_id].GridBound, self.exp_type[sim_id].GridBound))]
        else:
            [Time, GenSpaceTot] = self.read_gen_space_tot(sim_id, task_id, evol_pop, LD, time_type)
        ## LD = 1
        # make stable wild-type plot
        # colors of increasing shade [blue, yellow, green, orange], can later add: purple, grey from https://htmlcolorcodes.com/
        colors = [['#AED6F1', '#85C1E9', '#5DADE2', '#3498DB', '#2E86C1', '#2874A6', '#21618C', '#1B4F72'],
                  ['#F9E79F', '#F7DC6F', '#F4D03F', '#F1C40F', '#D4AC0D', '#B7950B', '#9A7D0A', '#7D6608'],
                  ['#ABEBC6', '#82E0AA', '#58D68D', '#2ECC71', '#28B463', '#239B56', '#1D8348', '#186A3B'],
                  ['#F5CBA7', '#F0B27A', '#EB984E', '#E67E22', '#CA6F1E', '#AF601A', '#935116', '#784212']]
        ## colors = [['#3498DB', '#2E86C1', '#2874A6', '#21618C', '#1B4F72', '#AED6F1', '#85C1E9', '#5DADE2'],
        ##          ['#F1C40F', '#D4AC0D', '#B7950B', '#9A7D0A', '#7D6608', '#AED6F1', '#85C1E9', '#5DADE2']]
        ax.tick_params(axis='y', labelcolor='tab:blue')
        ax.set_xticks([1, LD + 1, self.exp_type[sim_id].GridBound])
        ax.set_xticklabels(["$1$", "$R+1$", "$L$"], fontsize=12)
        ax.set_ylim([0, self.exp_type[sim_id].CarryingCapacity])
        ax.set_xlim([1, self.exp_type[sim_id].GridBound])

        # plot profile if a single swarming population is present
        if self.exp_type[sim_id].PopulationNumber == 1 and self.exp_type[sim_id].ConsiderSwarm[0]:
            if isinstance(self.exp_type[sim_id], SwarmDeadlyA):
                swarm_density = params[1]
                self.exp_type[sim_id].DeathRate[0] = params[2]
            elif isinstance(self.exp_type[sim_id], DensSwitchA):
                swarm_density = params[1]
                self.exp_type[sim_id].DeathRate[0] = params[0]
                if params[2] > params[3]:
                    colors[0], colors[1]=colors[1], colors[0]
            S = self.S(GenSpaceTot, 0, swarm_density, sim_id)
            for g in range(self.exp_type[sim_id].GenBound):
                g = self.exp_type[sim_id].GenBound - 1 - g
                N_sim = np.zeros(2*self.exp_type[sim_id].GridBound-1)
                for pos in range(self.exp_type[sim_id].GridBound):
                    for g_i in range(g + 1):
                        N_sim[2*pos] += GenSpaceTot[0][g_i][pos]
                    if pos > 0:
                        N_sim[2*pos-1] = (N_sim[2*pos]+N_sim[2*(pos-1)])/2
                if g == self.exp_type[sim_id].GenBound - 1:
                    if S == 0:
                        ax.fill_between([k/2 + 1 for k in range(2*self.exp_type[sim_id].GridBound-1)], N_sim, color=colors[0][g], label='')
                    elif S == self.exp_type[sim_id].GridBound:
                        ax.fill_between([k/2 + 1 for k in range(2*self.exp_type[sim_id].GridBound-1)], N_sim, color=colors[1][g], label='')
                    else:
                        ax.fill_between([k/2 + 1 for k in range(2*S)], N_sim[:2*S], color=colors[1][g], label='')
                        ax.fill_between([k/2 + S + 1/2 for k in range(2*self.exp_type[sim_id].GridBound-2*S)], N_sim[2*S-1:], color=colors[0][g], label='simulation')
                    if not analytics:
                        ax.plot([k/2 + 1 for k in range(2*self.exp_type[sim_id].GridBound-1)], N_sim, color='black')
                else:
                    if S == 0:
                        ax.fill_between([k/2 + 1 for k in range(2*self.exp_type[sim_id].GridBound-1)], N_sim, color=colors[0][g], label='')
                    elif S == self.exp_type[sim_id].GridBound:
                        ax.fill_between([k/2 + 1 for k in range(2*self.exp_type[sim_id].GridBound-1)], N_sim, color=colors[1][g], label='')
                    else:
                        ax.fill_between([k/2 + 1 for k in range(2*S)], N_sim[:2*S], color=colors[1][g], label='')
                        ax.fill_between([k/2 + S + 1/2 for k in range(2*self.exp_type[sim_id].GridBound-2*S)], N_sim[2*S-1:], color=colors[0][g], label='')
        # plot a profile otherwise
        else:
            for plot_pop in range(self.exp_type[sim_id].PopulationNumber):
                color = plot_pop
                plot_pop = self.exp_type[sim_id].PopulationNumber - 1 - plot_pop
                for g in range(self.exp_type[sim_id].GenBound):
                    g = self.exp_type[sim_id].GenBound - 1 - g
                    N_sim = np.zeros(self.exp_type[sim_id].GridBound)
                    for x in range(self.exp_type[sim_id].GridBound):
                        for g_i in range(g + 1):
                            for plot_pop_i in range(plot_pop + 1):
                                N_sim[x] += GenSpaceTot[plot_pop_i][g_i][x]
                    if g == self.exp_type[sim_id].GenBound - 1 and plot_pop == self.exp_type[sim_id].PopulationNumber - 1:
                        ax.fill_between([k + 1 for k in range(self.exp_type[sim_id].GridBound)], N_sim, color=colors[color][g], label='simulation')
                        if not analytics:
                            ax.plot([k + 1 for k in range(self.exp_type[sim_id].GridBound)], N_sim, color='black')
                    elif g == self.exp_type[sim_id].GenBound - 1:
                        ax.fill_between([k + 1 for k in range(self.exp_type[sim_id].GridBound)], N_sim, color=colors[color][g])
                        if not analytics:
                            ax.plot([k + 1 for k in range(self.exp_type[sim_id].GridBound)],  N_sim, linestyle='--', color='black')
                    else:
                        ax.fill_between([k + 1 for k in range(self.exp_type[sim_id].GridBound)], N_sim, color=colors[color][g], label='')

        # prepare data for r_net plot
        if show_r_net:
            plot_pop = 0    # assumes identical birth and death rates for all populations
            r_net = np.zeros(self.exp_type[sim_id].GridBound)
            N_tot = np.zeros(self.exp_type[sim_id].GridBound)
            for pos in range(self.exp_type[sim_id].GridBound):
                for pop in range(self.exp_type[sim_id].PopulationNumber):
                    for gen in range(self.exp_type[sim_id].GenBound):
                        N_tot[pos] += GenSpaceTot[pop][gen][pos]
                if pos < LD:
                    r_net[pos] = self.exp_type[sim_id].BirthRate[plot_pop]*(1-N_tot[pos]/self.exp_type[sim_id].CarryingCapacity)-self.exp_type[sim_id].DeathRate[plot_pop]
                elif pos < LD+1:
                    r_net[pos] = (self.exp_type[sim_id].BirthRate[plot_pop]+self.exp_type[sim_id].StressBirthRate[plot_pop]) * (1 - N_tot[pos] / self.exp_type[sim_id].CarryingCapacity) - (self.exp_type[sim_id].DeathRate[plot_pop]+self.exp_type[sim_id].StressDeathRate[plot_pop])
                else:
                    r_net[pos] = - (self.exp_type[sim_id].DeathRate[plot_pop]+self.exp_type[sim_id].StressDeathRate[plot_pop])
            # make the plot
            ax2 = ax.twinx()
            ax2.bar([k + 1 for k in range(self.exp_type[sim_id].GridBound)], r_net, color='tab:red')
            ax2.tick_params(axis='y', labelcolor='tab:red')
            ax2.set_ylim([-1.5 * self.exp_type[sim_id].DeathRate[plot_pop], self.exp_type[sim_id].BirthRate[plot_pop]-self.exp_type[sim_id].DeathRate[plot_pop]/2])
            ax2.set_yticks([])
            # make labels
            if r_labels:
                ax2.set_ylabel('mutant net birth rate', color='tab:red', fontsize=12)
                ax2.tick_params(axis='y', labelcolor='tab:red')
                yticks = [-self.exp_type[sim_id].DeathRate[plot_pop], 0, self.exp_type[sim_id].DeathRate[plot_pop]]
                ylabels = ['$-\delta$', '$0$', '$\delta$']
                for i in range(int((self.exp_type[sim_id].BirthRate[plot_pop]-self.exp_type[sim_id].DeathRate[plot_pop]/2)/self.exp_type[sim_id].DeathRate[plot_pop])):
                    if i > 0:
                        yticks.append((i + 1) * self.exp_type[sim_id].DeathRate[plot_pop])
                        ylabels.append(str(i + 1) + "$\\delta$")
                ax2.set_yticks(yticks)
                ax2.set_yticklabels(ylabels, fontsize=12)
        # make analytics
        if analytics:
            staircase = self.exp_type[sim_id].make_staircase(params, LD)
            N_plot = np.zeros(self.exp_type[sim_id].GridBound)
            for plot_pop in range(self.exp_type[sim_id].PopulationNumber):
                for pos in range(self.exp_type[sim_id].GridBound):
                    N_plot[pos] += staircase.N[plot_pop*self.exp_type[sim_id].GridBound+pos]
                if plot_pop < self.exp_type[sim_id].PopulationNumber-1:
                    ax.plot([i+1 for i in range(self.exp_type[sim_id].GridBound)], N_plot, linestyle='--', color='black', label='analytics')
                else:
                    ax.plot([i + 1 for i in range(self.exp_type[sim_id].GridBound)], N_plot, color='black', label='analytics')
        # make labels
        if pop_labels:
            K = self.exp_type[sim_id].CarryingCapacity
            ax.set_ylabel('$N_x$=$\\#$ cells', color='tab:blue', fontsize=12)
            ax.set_yticks([0, K / 4, K / 2, 3 * K / 4, K])
            ax.set_yticklabels(["0", "K/4", "K/2", "3K/4", "K"], fontsize=12)
        else:
            ax.set_yticks([])
        if x_labels:
            ax.set_xlabel('space x', fontsize=12)


###
# EXPERIMENT DATA TYPE
###


class MotSwitchA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 2
    # save population data
    InitCellsNum = [1e2, 1e2]
    InitCellsGen = [0, 0]
    InitCellsPos = [0, 0]
    BirthRate = [1, 1]
    StressBirthRate = [0, 0]
    ResistCost = [0, 0]
    DeathRate = [1e-1, 1e-1]
    StressDeathRate = [0, 0]
    MuteUp = [1e-7, 1e-7]
    MuteDown = [1e-4, 1e-4]
    StressMuteUp = [0, 0]
    StressMuteDown = [0, 0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5, 0.5]
    StressDepChemotax = [False, False]
    Chemokin = [0, 0]
    switch_params = [0, 1e-3, 1e-2, 1e-1, 1]
    switch_names = ["0", "1E3", "1E2", "1E1", "1"]
    ConsiderSwarm = [False, False]
    SwarmingDens = [5e4, 5e4]
    SwarmingMove = [1, 1]
    ConsiderHGT = [False, False]
    HGTRate = [0, 0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        switch_index = task_id // (self.repeat * mov_options)
        x = task_id % mov_options
        move0_index = 0
        while x >= move0_index * (move0_index + 1) / 2:
            move0_index += 1
        move0_index -= 1
        move1_index = x - move0_index * (move0_index + 1)/2
        switch_index = int(switch_index)
        move0_index = int(move0_index)
        move1_index = int(move1_index)
        id = "M" + self.move_names[move0_index] + "_M" + self.move_names[move1_index] + "_S" + self.switch_names[switch_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        move0_index = self.move_params.index(params[0])
        move1_index = self.move_params.index(params[1])
        switch_index = self.switch_params.index(params[2])
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        task_id = switch_index*self.repeat*mov_options+run_index*mov_options+move0_index*(move0_index+1)/2+move1_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        switch_index = task_id // (self.repeat * mov_options)
        x = task_id % mov_options
        move0_index = 0
        while x >= move0_index * (move0_index + 1) / 2:
            move0_index += 1
        move0_index -= 1
        move1_index = x - move0_index * (move0_index + 1)/2
        switch_index = int(switch_index)
        move0_index = int(move0_index)
        move1_index = int(move1_index)
        params = (self.move_params[move0_index], self.move_params[move1_index], self.switch_params[switch_index])
        run_index = int(int(task_id / mov_options) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp(), sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[1].Move = params[1]
        p_prop[0].SwitchDown = params[2]
        p_prop[1].SwitchDown = params[2]
        p_prop[0].SwitchUp = params[2]
        p_prop[1].SwitchUp = params[2]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        l_prop.PopulationNumber = 2
        return sc.Staircase(l_prop,p_prop)


class ChemotaxA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    chemotax_params = [1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1]
    chemotax_names = ["1E1", "2E1", "3E1", "4E1", "5E1", "6E1", "7E1", "8E1", "9E1"]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        chemotax_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_C" + self.chemotax_names[chemotax_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        chemotax_index = self.chemotax_params.index(params[1])
        task_id = chemotax_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        chemotax_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.chemotax_params[chemotax_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].Chemotax = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class StandardA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    death_params = [1e-1, 3e-1]
    death_names = ["1E1", "3E1"]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        death_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_D" + self.death_names[death_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        death_index = self.death_params.index(params[1])
        task_id = death_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        death_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.death_params[death_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].DeathRate = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class SpecialA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    special_params = [1e-2, 1e-2, 1e-1]
    special_names = ["RC1E2", "SB1E2", "SD1E1"]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        special_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_" + self.special_names[special_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        task_id = params[1] * self.repeat * mov_num + run_index * mov_num + move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        special_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], special_index)
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        if params[1] == 0:
            p_prop[0].ResistCost = self.special_params[0]
        elif params[1] == 1:
            p_prop[0].StressBirthRate = self.special_params[1]
        elif params[1] == 2:
            p_prop[0].StressDeathRate = self.special_params[2]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class ResistCostA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    cost_params = [1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 1.78e-1, 3.16e-1, 5.62e-1]
    cost_names = ["1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "1d78E1", "3E1", "5E1"]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-2, 1]
    move_names = ["1E2", "1"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        cost_num = len(self.cost_params)
        move_index = int(task_id / (self.repeat * cost_num))
        cost_index = int(task_id % cost_num)
        id = "M" + self.move_names[move_index] + "_RC" + self.cost_names[cost_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        cost_num = len(self.cost_params)
        cost_index = self.cost_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = move_index*self.repeat*cost_num+run_index*cost_num+cost_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        cost_num = len(self.cost_params)
        move_index = int(task_id / (self.repeat * cost_num))
        cost_index = int(task_id % cost_num)
        params = (self.move_params[move_index], self.cost_params[cost_index])
        run_index = int(int(task_id / cost_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].ResistCost = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class StressBirthA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    sbirth_params = [1e-3, 1e-2, 3.16e-2, 5.62e-2, 1e-1, 1.26e-1, 1.58e-1, 2e-1, 2.51e-1, 3.16e-1, 4e-1, 5.01e-1, 6.31e-1, 7.94e-1]
    sbirth_names = ["1E3", "1E2", "3E2", "5E2", "1E1", "1d26E1", "1d58E1", "2E1", "2d51E1", "3E1", "4E1", "5E1", "6E1", "7E1"]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-2, 1]
    move_names = ["1E2", "1"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        sbirth_num = len(self.sbirth_params)
        move_index = int(task_id / (self.repeat * sbirth_num))
        sbirth_index = int(task_id % sbirth_num)
        id = "M" + self.move_names[move_index] + "_SB" + self.sbirth_names[sbirth_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        sbirth_num = len(self.sbirth_params)
        sbirth_index = self.sbirth_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = move_index * self.repeat * sbirth_num + run_index * sbirth_num + sbirth_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        sbirth_num = len(self.sbirth_params)
        move_index = int(task_id / (self.repeat * sbirth_num))
        sbirth_index = int(task_id % sbirth_num)
        params = (self.move_params[move_index], self.sbirth_params[sbirth_index])
        run_index = int(int(task_id / sbirth_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].StressBirthRate = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class StressDeathA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    sdeath_params = [0, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 1]
    sdeath_names = ["0", "1E1", "2E1", "3E1", "4E1", "5E1", "6E1", "7E1", "8E1", "9E1", "1"]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-2, 1]
    move_names = ["1E2", "1"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        sdeath_num = len(self.sdeath_params)
        move_index = int(task_id / (self.repeat * sdeath_num))
        sdeath_index = int(task_id % sdeath_num)
        id = "M" + self.move_names[move_index] + "_SD" + self.sdeath_names[sdeath_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        sdeath_num = len(self.sdeath_params)
        sdeath_index = self.sdeath_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = move_index * self.repeat * sdeath_num + run_index * sdeath_num + sdeath_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        sdeath_num = len(self.sdeath_params)
        move_index = int(task_id / (self.repeat * sdeath_num))
        sdeath_index = int(task_id % sdeath_num)
        params = (self.move_params[move_index], self.sdeath_params[sdeath_index])
        run_index = int(int(task_id / sdeath_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].StressDeathRate = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class SpecialB:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    special_params = [1e-2, 1e-2, 10, 1] #(1) resitance cost, (2) stress birth rate, (3) stress death rate + stress birth rate
    special_names = ["RC1E2", "SB1E2", "SB1_SD10"]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        special_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_" + self.special_names[special_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        task_id = params[1] * self.repeat * mov_num + run_index * mov_num + move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        special_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], special_index)
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        if params[1] == 0:
            p_prop[0].ResistCost = self.special_params[0]
        elif params[1] == 1:
            p_prop[0].StressBirthRate = self.special_params[1]
        elif params[1] == 2:
            p_prop[0].StressDeathRate = self.special_params[2]
            p_prop[0].StressBirthRate = self.special_params[3]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class ResistCostB:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    cost_params = [1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 1.78e-1, 3.16e-1, 5.62e-1, 6.3e-1, 7.94e-1]
    cost_names = ["1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "1d78E1", "3E1", "5E1", "6E1", "7E1"]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-2, 1]
    move_names = ["1E2", "1"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        cost_num = len(self.cost_params)
        move_index = int(task_id / (self.repeat * cost_num))
        cost_index = int(task_id % cost_num)
        id = "M" + self.move_names[move_index] + "_RC" + self.cost_names[cost_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        cost_num = len(self.cost_params)
        cost_index = self.cost_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = move_index*self.repeat*cost_num+run_index*cost_num+cost_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        cost_num = len(self.cost_params)
        move_index = int(task_id / (self.repeat * cost_num))
        cost_index = int(task_id % cost_num)
        params = (self.move_params[move_index], self.cost_params[cost_index])
        run_index = int(int(task_id / cost_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].ResistCost = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class StressBirthDeathA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [1]
    ResistCost = [0]
    DeathRate = [1e-1]
    sdeath_params = [0.9, 1.16, 1.48, 1.9, 2.41, 3.06, 3.88, 4.91, 6.2, 7.84, 9.9]
    sdeath_names = ["9E1", "1d16", "1d48", "1d9", "2d41", "3d06", "3d88", "4", "6", "7", "9"]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-2, 1]
    move_names = ["1E2", "1"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        sdeath_num = len(self.sdeath_params)
        move_index = int(task_id / (self.repeat * sdeath_num))
        sdeath_index = int(task_id % sdeath_num)
        id = "M" + self.move_names[move_index] + "_SB1_SD" + self.sdeath_names[sdeath_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        sdeath_num = len(self.sdeath_params)
        sdeath_index = self.sdeath_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = move_index * self.repeat * sdeath_num + run_index * sdeath_num + sdeath_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        sdeath_num = len(self.sdeath_params)
        move_index = int(task_id / (self.repeat * sdeath_num))
        sdeath_index = int(task_id % sdeath_num)
        params = (self.move_params[move_index], self.sdeath_params[sdeath_index])
        run_index = int(int(task_id / sdeath_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].StressDeathRate = params[1]
        p_prop[0].StressBirthRate = self.StressBirthRate[0]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class StressDeathB:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    sdeath_params = [0, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1]
    sdeath_names = ["0", "1E1", "2E1", "3E1", "4E1", "5E1", "6E1", "7E1", "8E1"]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        sdeath_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_SD" + self.sdeath_names[sdeath_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        sdeath_index = self.sdeath_params.index(params[1])
        move_index = self.move_params.index(params[0])
        task_id = sdeath_index * self.repeat * mov_num + run_index * mov_num + move_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        sdeath_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.sdeath_params[sdeath_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # make staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].StressDeathRate = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class SwitchA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 2
    # save population data
    InitCellsNum = [1e2, 1e2]
    InitCellsGen = [0, 0]
    InitCellsPos = [0, 0]
    BirthRate = [1, 1]
    StressBirthRate = [0, 0]
    ResistCost = [0, 0]
    death_params = [1e-1, 3e-1]
    death_names = ["1E1", "3E1"]
    StressDeathRate = [0, 0]
    MuteUp = [1e-7, 1e-7]
    MuteDown = [1e-4, 1e-4]
    StressMuteUp = [0, 0]
    StressMuteDown = [0, 0]
    Move = [1e-3, 10]
    Chemotax = [0.5, 0.5]
    StressDepChemotax = [False, False]
    Chemokin = [0, 0]
    switch_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    switch_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    ConsiderSwarm = [False, False]
    SwarmingDensity = [5e4, 5e4]
    SwarmingMove = [1, 1]
    ConsiderHGT = [False, False]
    HGTRate = [0, 0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        switch_num = len(self.switch_params)
        death_index = int(task_id / (self.repeat * switch_num))
        switch_index = int(task_id % switch_num)
        id = "S" + self.switch_names[switch_index] + "_D" + self.death_names[death_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        switch_num = len(self.switch_params)
        switch_index = self.switch_params.index(params[0])
        death_index = self.death_params.index(params[1])
        task_id = death_index * self.repeat * switch_num + run_index * switch_num + switch_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        switch_num = len(self.switch_params)
        death_index = int(task_id / (self.repeat * switch_num))
        switch_index = int(task_id % switch_num)
        params = (self.switch_params[switch_index], self.death_params[death_index])
        run_index = int(int(task_id / switch_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp(), sc.PopProp()]
        p_prop[0].SwitchDown = params[0]
        p_prop[1].SwitchDown = params[0]
        p_prop[0].SwitchUp = params[0]
        p_prop[1].SwitchUp = params[0]
        p_prop[0].DeathRate = params[1]
        p_prop[0].DeathRate = params[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        l_prop.PopulationNumber = 2
        return sc.Staircase(l_prop,p_prop)


class MotSwitchDeadlyA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 2
    # save population data
    InitCellsNum = [1e2, 1e2]
    InitCellsGen = [0, 0]
    InitCellsPos = [0, 0]
    BirthRate = [1, 1]
    StressBirthRate = [0, 0]
    ResistCost = [0, 0]
    DeathRate = [3e-1, 3e-1]
    StressDeathRate = [0, 0]
    MuteUp = [1e-7, 1e-7]
    MuteDown = [1e-4, 1e-4]
    StressMuteUp = [0, 0]
    StressMuteDown = [0, 0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5, 0.5]
    StressDepChemotax = [False, False]
    Chemokin = [0, 0]
    switch_params = [0, 1e-3, 1e-1, 4e-1, 5]
    switch_names = ["0", "1E3", "1E1", "4E1", "5"]
    ConsiderSwarm = [False, False]
    SwarmingDensity = [5e4, 5e4]
    SwarmingMove = [1, 1]
    ConsiderHGT = [False, False]
    HGTRate = [0, 0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        switch_index = task_id // (self.repeat * mov_options)
        x = task_id % mov_options
        move0_index = 0
        while x >= move0_index * (move0_index + 1) / 2:
            move0_index += 1
        move0_index -= 1
        move1_index = x - move0_index * (move0_index + 1)/2
        switch_index = int(switch_index)
        move0_index = int(move0_index)
        move1_index = int(move1_index)
        id = "M" + self.move_names[move0_index] + "_M" + self.move_names[move1_index] + "_S" + self.switch_names[switch_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        move0_index = self.move_params.index(params[0])
        move1_index = self.move_params.index(params[1])
        switch_index = self.switch_params.index(params[2])
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        task_id = switch_index*self.repeat*mov_options+run_index*mov_options+move0_index*(move0_index+1)/2+move1_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        mov_options = mov_num * (mov_num + 1) / 2
        switch_index = task_id // (self.repeat * mov_options)
        x = task_id % mov_options
        move0_index = 0
        while x >= move0_index * (move0_index + 1) / 2:
            move0_index += 1
        move0_index -= 1
        move1_index = x - move0_index * (move0_index + 1)/2
        switch_index = int(switch_index)
        move0_index = int(move0_index)
        move1_index = int(move1_index)
        params = (self.move_params[move0_index], self.move_params[move1_index], self.switch_params[switch_index])
        run_index = int(int(task_id / mov_options) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp(), sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[1].Move = params[1]
        p_prop[0].SwitchDown = params[2]
        p_prop[1].SwitchDown = params[2]
        p_prop[0].SwitchUp = params[2]
        p_prop[1].SwitchUp = params[2]
        p_prop[0].DeathRate = self.DeathRate[0]
        p_prop[1].DeathRate = self.DeathRate[1]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        l_prop.PopulationNumber = 2
        return sc.Staircase(l_prop, p_prop)


class SwarmLagA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    init_cell_params = [1e2, 2.15e2, 4.64e2, 1e3, 2.15e3, 4.64e3, 1e4, 2.15e4, 4.64e4]
    init_cell_names = ["1E2", "2E2", "4E2", "1E3", "2E3", "4E3", "1E4", "2E4", "4E4"]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    Move = [1e-3]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [True]
    swarm_move_params = [3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    swarm_move_names = ["3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    SwarmingDensity = [5e4]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        move_num = len(self.swarm_move_params)
        init_cell_index = int(task_id / (self.repeat * move_num))
        move_index = int(int(task_id) % move_num)
        id="SwM" + self.swarm_move_names[move_index] + "_IC" + self.init_cell_names[init_cell_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        move_num = len(self.swarm_move_params)
        move_index = self.swarm_move_params.index(params[0])
        init_cell_index = self.init_cell_params.index(params[1])
        task_id = init_cell_index*self.repeat*move_num + run_index*move_num + move_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        move_num = len(self.swarm_move_params)
        init_cell_index = int(task_id / (self.repeat * move_num))
        move_index = int(int(task_id) % move_num)
        params = (self.swarm_move_params[move_index], self.init_cell_params[init_cell_index])
        run_index = int(int(task_id/move_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].SwarmingMove = params[0]
        p_prop[0].InitCellsNum = params[1]
        p_prop[0].ConsiderSwarm = True
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class SwarmDeadlyA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    death_params = [1e-1, 3e-1]
    death_names = ["1E1", "3E1"]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    Move = [1e-3]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [True]
    swarm_density_params = [1E2, 2.15E2, 4.64E2, 1E3, 2.15E3, 4.64E3, 1E4, 2.15E4, 4.64E4, 6E4, 7E4, 8E4, 9E4]
    swarm_density_names = ["1E2", "2E2", "4E2", "1E3", "2E3", "4E3", "1E4", "2E4", "4E4", "6E4", "7E4", "8E4", "9E4"]
    swarm_move_params = [3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10]
    swarm_move_names = ["3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        dense_num = len(self.swarm_density_params)
        move_num = len(self.swarm_move_params)
        death_index = int(task_id / (self.repeat * move_num * dense_num))
        dense_index = int(int(task_id % (move_num * dense_num)) / move_num)
        move_index = int(int(task_id % (move_num * dense_num)) % move_num)
        id = "SwM"+self.swarm_move_names[move_index]+"_SwD"+self.swarm_density_names[dense_index]+"_D"+self.death_names[death_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        dense_num = len(self.swarm_density_params)
        move_num = len(self.swarm_move_params)
        move_index = self.swarm_move_params.index(params[0])
        density_index = self.swarm_density_params.index(params[1])
        death_index = self.death_params.index(params[2])
        task_id = death_index*self.repeat*dense_num*move_num + run_index*dense_num*move_num + density_index*move_num + move_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        dense_num = len(self.swarm_density_params)
        move_num = len(self.swarm_move_params)
        death_index = int(task_id / (self.repeat * move_num * dense_num))
        dense_index = int(int(task_id % (move_num * dense_num)) / move_num)
        move_index = int(int(task_id % (move_num * dense_num)) % move_num)
        params = (self.swarm_move_params[move_index], self.swarm_density_params[dense_index], self.death_params[death_index])
        run_index = int(int(task_id / (move_num * dense_num)) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].SwarmingMove = params[0]
        p_prop[0].SwarmingDensity = params[1]
        p_prop[0].DeathRate = params[2]
        p_prop[0].ConsiderSwarm = True
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class HGTA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [True]
    hgt_params = [1e-7, 3.16e-7, 1e-6, 3.16e-6, 1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3]
    hgt_names = ["1E7", "3E7", "1E6", "3E6", "1E5", "3E5", "1E4", "3E4", "1E3"]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        hgt_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_HGT" + self.hgt_names[hgt_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        hgt_index = self.hgt_params.index(params[1])
        task_id = hgt_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        hgt_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.hgt_params[hgt_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].HGTRate = params[1]
        p_prop[0].ConsiderHGT = True
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class HGTDeadlyA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    ResistCost = [0]
    DeathRate = [3e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [True]
    hgt_params = [1e-7, 3.16e-7, 1e-6, 3.16e-6, 1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3]
    hgt_names = ["1E7", "3E7", "1E6", "3E6", "1E5", "3E5", "1E4", "3E4", "1E3"]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        hgt_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_HGT" + self.hgt_names[hgt_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        hgt_index = self.hgt_params.index(params[1])
        task_id = hgt_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        hgt_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.hgt_params[hgt_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].HGTRate = params[1]
        p_prop[0].ConsiderHGT = True
        p_prop[0].DeathRate = self.DeathRate[0]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class DensSwitchA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [50]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    DeathRate = [None]
    StressBirthRate = [0]
    ResistCost = [0]
    death_params = [1e-1, 3e-1]
    death_names = ["1E1", "3E1"]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1E-5, 3.16E-5, 1E-4, 3.16E-4, 1E-3, 3.16E-3, 1E-2, 3.16E-2, 1E-1, 3.16E-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [True]
    switch_density_params = [1E2, 1E3, 1E4, 3.14E4, 9.5E4]
    switch_density_names = ["1E2", "1E3", "1E4", "3E4", "9E4"]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        dense_num = len(self.switch_density_params)
        move_num = len(self.move_params)
        death_index = int(task_id / (self.repeat * move_num * move_num * dense_num))
        dense_index = int(int(task_id % (move_num * move_num * dense_num)) / (move_num*move_num))
        move_low_dens_index = int(int(int(task_id % (move_num * move_num * dense_num)) % (move_num * move_num)) / move_num)
        move_high_dens_index = int(int(int(task_id % (move_num * move_num * dense_num)) % (move_num * move_num)) % move_num)
        id = "M"+self.move_names[move_low_dens_index]+"_SwM"+self.move_names[move_high_dens_index]+"_SwD"+self.switch_density_names[dense_index]+"_D"+self.death_names[death_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        dense_num = len(self.switch_density_params)
        move_num = len(self.move_params)
        death_index = self.death_params.index(params[0])
        dense_index = self.switch_density_params.index(params[1])
        move_low_dens_index = self.move_params.index(params[2])
        move_high_dens_index = self.move_params.index(params[3])
        task_id = death_index*self.repeat*dense_num*move_num*move_num + run_index*dense_num*move_num*move_num + dense_index*move_num*move_num + move_low_dens_index*move_num + move_high_dens_index
        return int(task_id)

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        dense_num = len(self.switch_density_params)
        move_num = len(self.move_params)
        death_index = int(task_id / (self.repeat * move_num * move_num * dense_num))
        dense_index = int(int(task_id % (move_num * move_num * dense_num)) / (move_num*move_num))
        move_low_dens_index = int(int(int(task_id % (move_num * move_num * dense_num)) % (move_num * move_num)) / move_num)
        move_high_dens_index = int(int(int(task_id % (move_num * move_num * dense_num)) % (move_num * move_num)) % move_num)
        params = (self.death_params[death_index], self.switch_density_params[dense_index], self.move_params[move_low_dens_index], self.move_params[move_high_dens_index])
        run_index = int(int(task_id / (move_num * move_num * dense_num)) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].InitCellsNum = 50
        p_prop[0].DeathRate = params[0]
        p_prop[0].SwarmingDensity = params[1]
        p_prop[0].Move = params[2]
        p_prop[0].SwarmingMove = params[3]
        p_prop[0].ConsiderSwarm = True
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class ResCostA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    cost_params = [3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1]
    cost_names = ["3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1"]
    DeathRate = [1e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        cost_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_RC" + self.cost_names[cost_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        cost_index = self.cost_params.index(params[1])
        task_id = cost_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        cost_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.cost_params[cost_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].ResistCost = params[1]
        p_prop[0].DeathRate = self.DeathRate[0]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)


class ResCostDeadlyA:
    repeat = 10
    adap_num = 7
    # save lattice data
    D = 1
    CarryingCapacity = 1e5
    GenBound = 8
    GridBound = 8
    PopulationNumber = 1
    # save population data
    InitCellsNum = [1e2]
    InitCellsGen = [0]
    InitCellsPos = [0]
    BirthRate = [1]
    StressBirthRate = [0]
    cost_params = [3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1]
    cost_names = ["3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1"]
    DeathRate = [3e-1]
    StressDeathRate = [0]
    MuteUp = [1e-7]
    MuteDown = [1e-4]
    StressMuteUp = [0]
    StressMuteDown = [0]
    move_params = [1e-5, 3.16e-5, 1e-4, 3.16e-4, 1e-3, 3.16e-3, 1e-2, 3.16e-2, 1e-1, 3.16e-1, 1, 3.16, 10]
    move_names = ["1E5", "3E5", "1E4", "3E4", "1E3", "3E3", "1E2", "3E2", "1E1", "3E1", "1", "3", "10"]
    Chemotax = [0.5]
    StressDepChemotax = [False]
    Chemokin = [0]
    SwitchUp = [0]
    SwitchDown = [0]
    ConsiderSwarm = [False]
    SwarmingDensity = [5e4]
    SwarmingMove = [1]
    ConsiderHGT = [False]
    HGTRate = [0]

    # Find name_ID from task_ID
    def taskIDtoID(self, task_id):
        mov_num = len(self.move_params)
        cost_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        id = "M" + self.move_names[move_index] + "_RC" + self.cost_names[cost_index]
        return id

    # Find taskID from params and run_index(=0,...,repeat)
    def IDtoTaskID(self, params, run_index):
        mov_num = len(self.move_params)
        move_index = self.move_params.index(params[0])
        cost_index = self.cost_params.index(params[1])
        task_id = cost_index*self.repeat*mov_num+run_index*mov_num+move_index
        return task_id

    # Find params from task_ID
    def taskIDtoParams(self, task_id):
        mov_num = len(self.move_params)
        cost_index = int(task_id / (self.repeat * mov_num))
        move_index = int(task_id % mov_num)
        params = (self.move_params[move_index], self.cost_params[cost_index])
        run_index = int(int(task_id / mov_num) % self.repeat)
        return [params, run_index]

    # Create staircase object
    def make_staircase(self, params, LD):
        p_prop = [sc.PopProp()]
        p_prop[0].Move = params[0]
        p_prop[0].ResistCost = params[1]
        p_prop[0].DeathRate = self.DeathRate[0]
        l_prop = sc.LatProp()
        l_prop.LD = LD
        return sc.Staircase(l_prop, p_prop)
