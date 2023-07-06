import random
import dimod
import dwave.inspector
import csv
import math
import dwave_networkx as dnx
import networkx as nx
import time
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite, LazyFixedEmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
from dimod.vartypes import SPIN
from dimod import SampleSet
from minorminer import find_embedding
from dwave.embedding.chain_strength import uniform_torque_compensation as UTC
from dwave.embedding import embed_bqm, unembed_sampleset, EmbeddedStructure
from functools import partial
from copy import deepcopy
from dwave.samplers import SimulatedAnnealingSampler


SASampler = SimulatedAnnealingSampler()
# CSampler = DWaveSampler(solver={'topology__type': 'pegasus'})

from helpers import HillClimbChimeraSampler


def fe(S, T, **kwargs):
    print('find_embedding called')
    return find_embedding(S, T, random_seed=123123)

FCSampler = HillClimbChimeraSampler()

# C_Composite = LazyFixedEmbeddingComposite(CSampler, find_embedding=fe)
FC_Composite = LazyFixedEmbeddingComposite(FCSampler, find_embedding=fe)





print('reading input')

iname: str = input()

n = int(input())

model = BinaryQuadraticModel(vartype=SPIN)

for i in range(n):
    model.set_linear(i,int(input()))

while(True):
    line = (input())
    listt = line.split(' ')
    # print(line)
    # print(listt)
    temp = list(map(int,listt))
    if (temp == [-1]):
        break
    # print(temp)
    u,v,a = temp
    model.set_quadratic(u,v,a)

optimal_solution = int(input())


print(f'Optimal solution: {optimal_solution}')




EPS = 1e-9
cached : dict = {}


def compute(j: float, composite: dimod.Sampler, runs: int) -> dict[str,float]:
    print('Final run:')
    print(f'strength: {j}')
    
    print(time.time())
    sample_set = composite.sample(model, num_reads = runs, chain_strength = j)
    print(time.time())

    s1_cnt: float = 0
    s2_cnt: float = 0
    ncb_cnt: float = 0
    avg: float = 0
    best: float = 0
    for r in sample_set.record:           
        if ( r.energy < optimal_solution * 0.995 + EPS):
            s1_cnt += r.num_occurrences
                                
        if ( r.energy < optimal_solution * 0.97 + EPS):
            s2_cnt += r.num_occurrences
        avg += r.energy * r.num_occurrences
        best = min(best, r.energy)
        ncb_cnt += r.chain_break_fraction
    
    res: dict[str,float] = {}
    res['p5no'] = s1_cnt / runs
    res['3no'] = s2_cnt / runs
    res['avg'] = avg / runs
    res['break'] = ncb_cnt / runs
    res['best'] = best
    print(res)
    return res

DEF_D: float = 0.01
def get_d(j: float, delta: float) -> float:
    res = compute(j + delta / 2, FC_Composite, 1000)
    res_d = compute(j - delta / 2, FC_Composite, 1000)
    
    return math.log(res_d['break'] / res['break']) / delta


DEF_LB: float = 0.002
DEF_UB: float = 0.1
DEF_BS_RUNS: int = 15

def get_cs_range(lb: float, ub: float) -> tuple:
    start_low = 0.01
    start_high: float = 2 * UTC(model)
    res_lb : float = 0
    res_ub : float = 0
    l: float = 0
    r: float = 0
    
    l = start_low
    r = start_high
    for i in range(DEF_BS_RUNS):
        mid = (l+r)/2
        if (compute(mid, FC_Composite, 200)['break'] < lb):
            r = mid
        else:
            l = mid
    res_lb = (l+r)/2
    
    
    l = start_low
    r = start_high
    for i in range(DEF_BS_RUNS):
        mid = (l+r)/2
        if (compute(mid, FC_Composite, 200)['break'] < ub):
            r = mid
        else:
            l = mid
    res_ub = (l+r)/2
    print(f'result: {res_lb}, {res_ub}')
    return (res_lb, res_ub)

DEF_STEP: float = 0.002
DEF_TRIES: int = 12
DEF_SPLIT: int = 4
DEF_ROUNDING: float = 0.1

def hill_climb(lb: float, ub: float) -> list[float]:
    result: list[point] = []
    for tri in range(DEF_TRIES):
        point = random.random()*(ub-lb)/DEF_SPLIT+lb+(ub-lb)/DEF_SPLIT*(tri%DEF_SPLIT)
        # while True:
        #     point = random.random()*(ub-lb)+lb
        #     if (compute(point, FC_Composite, 200)['avg'] < optimal_solution * 0.65):
        #         break

        while True:
            stop = True
            u = DEF_STEP * (ub-lb)
            if get_d(point, DEF_D*(ub-lb)) < get_d(point + u, DEF_D*(ub-lb)):
                point += u
                stop = False
                break
            if get_d(point, DEF_D*(ub-lb)) < get_d(point - u, DEF_D*(ub-lb)):
                point -= u
                stop = False
                break
            if (stop):
                break
        stop = False
        for g in result:
            if (abs(g-point) <= DEF_ROUNDING):
                stop = True
        if (not stop):
            result.append(point)
    return result

start_time = time.time()

cs_range = get_cs_range(DEF_LB, DEF_UB)

result = hill_climb(cs_range[0], cs_range[1])
result.sort()

f = open(f'findgap2_runtime.txt', mode='a')
print(f'Run time: {(time.time() - start_time)} seconds')
f.write(f'{iname},{(time.time() - start_time)}\n')
f.close()

f = open(f'findgap2_results_p/{iname}.csv', 'w', encoding='utf-8')
f.write('abs,rel,p5no,3no,best,avg,break\n')
sus: float = UTC(model)
def output(j):
    print(j)
    print(j/sus)
    # u = compute(j, C_Composite, 1000)
    # print(u)
    # f.write(f'{j},{j/sus},{u["p5no"]},{u["3no"]},{u["best"]},{u["avg"]},{u["break"]}\n')

for cs in result:
    output(cs)
output(sus)

f.close()
print(UTC(model))   