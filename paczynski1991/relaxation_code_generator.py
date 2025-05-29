# This is a code which will use sympy to autogenerate a relaxation code for a fairly general class of
# relaxation method problems (including eigenproblems).
import sympy as sym
import os
import shutil

from sympy.codegen.ast import Assignment
import pdb

# compile and run with:
# /usr/bin/clang++ -Wall -std=c++11 main.cpp diffeq.cpp -o main; ./main

# Directory parameters
dirname = 'relaxation_paczynski1991/'
nr3dir = 'resources/'

outdir = f'{dirname}output/'

precision = 10
thresh = 1e-3

# Solver parameters
N_grid = 10

# Important symbol definitions
x = sym.Symbol('x', real=True)
x_min = 0.8
x_max = 10.

dx = sym.Symbol('dx', real=True)
xbar = sym.Symbol('xbar', real=True)

# Iteration parameters
n = sym.Symbol('n', real=True)
omega_spin = sym.Symbol('omega_spin', real=True)
a_inv = sym.Symbol('a_inv', real=True)

paramlist = [a_inv, omega_spin, n]
paramdict_type = {a_inv: 'Doub', omega_spin: 'Doub', n: 'Int'} # list: dtype, min, max, spacing
paramdict_min = {a_inv: 1e6, omega_spin: 0.8, n: 1}
paramdict_max = {a_inv: 2e6, omega_spin: 0.9, n: 1}
paramdict_points = {a_inv: 2, omega_spin: 2, n: None}

# Unknown functions to solve for
y = sym.Function('y', real=True)
omega = sym.Function('omega', real=True)
j_star = sym.Function('j_star', real=True)

Unknowns = [omega, y, j_star]

Unknowns_flat = {y: False, omega: False, j_star: True} # if flat, means just write out the first item (since they're all the same)

# Add initial guesses as strings # TODO
Guesses = ['pow(x_grid[ii], -1.5)',
           'pow(a_inv_ii / 1.5, 1./(2.*n+3.)) * pow(x_grid[ii], (6.*n+3.) / (4.*n+6.))',
           '0.'] # initial guesses

# Boundary conditions (as expressions which must vanish at the specified endpoint) - don't put derivatives into these
Left1 = omega_spin - omega(x)
Left2 = ((1 - 0.5 * omega_spin ** 2 * x_min ** 2) ** -2. - x_min ** 2) ** 0.5 - y(x)

Boundary_left = [Left1, Left2]

Right1 = x_max ** -1.5 - omega(x)

Boundary_right = [Right1]

# Equations, written as expressions which will be set to zero
Equation1 = x ** (3. * n - 1.5) * (j_star(x) - omega(x) * x ** 2) * a_inv / y(x) ** (2. * n + 3.) - sym.diff(omega(x), x)
Equation2 = omega(x) ** 2 * (x ** 2 + y(x) ** 2) ** 1.5 * x / y(x) - x / y(x) - sym.diff(y(x), x)
Equation3 = sym.diff(j_star(x), x) # const. j_star

Equations = [Equation1, Equation2, Equation3]

# "Bump parameters" - perturbing the initial guess each time
bumps = [0, 0, 0]

# Numerical Recipes parameters
run_itmax = 10
run_conv = 1.e-15
run_slowc = 1.
run_scalv = [1., 1., 1.]

#####################
### THE BLACK BOX ###
#####################

# Copy files
os.mkdir(dirname)
os.mkdir(outdir)

shutil.copy(f'{nr3dir}nr3.h', f'{dirname}nr3.h')
shutil.copy(f'{nr3dir}solvde.h', f'{dirname}solvde.h')

# Discretize the LEFT boundary conditions and find the relevant derivatives
s_matrix_assignments = []

for ii, Left in enumerate(Boundary_left):
    Equation_dummy = Left
    vars_endpoint = []
    
    for jj, param in enumerate(paramlist):
        var_symb = sym.Symbol(f'var_{param.name}')
        Equation_dummy = Equation_dummy.subs(param, var_symb)
    
    for jj, var in enumerate(Unknowns):
        var_endpoint = sym.Symbol(f'y[{jj}][0]')
        vars_endpoint.append(var_endpoint)
        
        Equation_dummy = Equation_dummy.subs(var(x), var_endpoint)
    
    Equation_dummy = Equation_dummy.subs(x, xbar)
    
    # Get discretized equations for s matrix elements
    first_index = len(Unknowns) - len(Boundary_left)
    
    for jj, var in enumerate(Unknowns):
        var_endpoint = vars_endpoint[jj]
        
        assignment_var = sym.Symbol(f's[{first_index+ii}][jsf]')
        assignment = Assignment(assignment_var, Equation_dummy)
        s_matrix_assignments.append(sym.ccode(assignment))
        
        assignment_var = sym.Symbol(f's[{first_index+ii}][indexv[{jj}]+{len(Unknowns)}]')
        assignment = Assignment(assignment_var, sym.diff(Equation_dummy, var_endpoint))
        s_matrix_assignments.append(sym.ccode(assignment))

left_assignments = s_matrix_assignments

# Discretize the RIGHT boundary conditions and find the relevant derivatives
s_matrix_assignments = []

for ii, Right in enumerate(Boundary_right):
    Equation_dummy = Right
    vars_endpoint = []
    
    for jj, param in enumerate(paramlist):
        var_symb = sym.Symbol(f'var_{param.name}')
        Equation_dummy = Equation_dummy.subs(param, var_symb)
    
    for jj, var in enumerate(Unknowns):
        var_endpoint = sym.Symbol(f'y[{jj}][mpt-1]')
        vars_endpoint.append(var_endpoint)

        Equation_dummy = Equation_dummy.subs(var(x), var_endpoint)
    
    Equation_dummy = Equation_dummy.subs(x, x_max)
    
    # Get discretized equations for s matrix elements    
    for jj, var in enumerate(Unknowns):
        var_endpoint = vars_endpoint[jj]
        
        assignment_var = sym.Symbol(f's[{ii}][jsf]')
        assignment = Assignment(assignment_var, Equation_dummy)
        s_matrix_assignments.append(sym.ccode(assignment))
        
        assignment_var = sym.Symbol(f's[{ii}][indexv[{jj}]+{len(Unknowns)}]')
        assignment = Assignment(assignment_var, sym.diff(Equation_dummy, var_endpoint))
        s_matrix_assignments.append(sym.ccode(assignment))

right_assignments = s_matrix_assignments

# Discretize the equations and find the relevant derivatives
s_matrix_assignments = []

for ii, Equation in enumerate(Equations):
    Equation_dummy = Equation
    
    vars_upper = []
    vars_lower = []
    
    for jj, param in enumerate(paramlist):
        var_symb = sym.Symbol(f'var_{param.name}')
        Equation_dummy = Equation_dummy.subs(param, var_symb)
    
    for jj, var in enumerate(Unknowns):
        var_upper = sym.Symbol(f'y[{jj}][k]')
        var_lower = sym.Symbol(f'y[{jj}][k-1]')
        
        vars_upper.append(var_upper)
        vars_lower.append(var_lower)
        
        Equation_dummy = Equation_dummy.subs(sym.diff(var(x), x), (var_upper - var_lower) / dx)
        Equation_dummy = Equation_dummy.subs(var(x), (var_upper + var_lower) / 2)
    
    Equation_dummy = Equation_dummy.subs(x, xbar)
    
    # Get discretized equations for s matrix elements
    for jj, var in enumerate(Unknowns):
        var_upper = vars_upper[jj]
        var_lower = vars_lower[jj]
        
        assignment_var = sym.Symbol(f's[{ii}][jsf]')
        assignment = Assignment(assignment_var, Equation_dummy)
        s_matrix_assignments.append(sym.ccode(assignment))
        
        assignment_var = sym.Symbol(f's[{ii}][indexv[{jj}]]')
        assignment = Assignment(assignment_var, sym.diff(Equation_dummy, var_lower))
        s_matrix_assignments.append(sym.ccode(assignment))
        
        assignment_var = sym.Symbol(f's[{ii}][indexv[{jj}]+{len(Unknowns)}]')
        assignment = Assignment(assignment_var, sym.diff(Equation_dummy, var_upper))
        s_matrix_assignments.append(sym.ccode(assignment))

interior_assignments = s_matrix_assignments

#############################
### WRITE DIFFEQ.HPP FILE ###
#############################
text = ''
text += '#ifndef diffeq_h\n'
text += '#define diffeq_h\n'
text += '\n'
text += 'struct Diffeq {\n'

for ii, key in enumerate(paramlist):
    text += f'    {paramdict_type[key]} var_{key};\n'    
text += '    VecDoub var_x_grid;\n'
text += '    Int mpt;\n'
text += '    Doub xbar, dx;\n'

text += '    Diffeq('
for ii, key in enumerate(paramlist):
    text += f'const {paramdict_type[key]} &in_{key}, '
text += 'const VecDoub in_x_grid);\n'
text += '    void smatrix(const Int k, const Int k1, const Int k2, const Int jsf, const Int is1, const Int isf, VecInt_I &indexv, MatDoub_O &s, MatDoub_I &y);\n'
text += '};\n'
text += '\n'
text += '#endif\n'

f = open(f'{dirname}diffeq.hpp', 'w')
f.write(text)
f.close()

#############################
### WRITE DIFFEQ.CPP FILE ###
#############################
text = ''
text += '#include <cmath>\n'
text += '#include <vector>\n'
text += '#include "nr3.h"\n'
text += '#include "diffeq.hpp"\n'
text += '\n'
text += 'using namespace std;\n'
text += '\n'

text += 'Diffeq::Diffeq('
for ii, param in enumerate(paramlist):
    text += f'const {paramdict_type[param]} &in_{param}, '
text += 'const VecDoub in_x_grid) {\n'

for ii, param in enumerate(paramlist):
    text += f'    var_{param} = in_{param};\n'

text += '    var_x_grid = in_x_grid;\n'
text += '    mpt = var_x_grid.size();\n'
text += '}\n'
text += '\n'

text += 'void Diffeq::smatrix(const Int k, const Int k1, const Int k2, const Int jsf, const Int is1, const Int isf, VecInt_I &indexv, MatDoub_O &s, MatDoub_I &y) {\n'
text += '    dx = var_x_grid[k] - var_x_grid[k-1];\n'
text += '    xbar = 0.5 * (var_x_grid[k] + var_x_grid[k-1]);\n'
text += '    \n'

text += '    if (k == k1) {\n'
for ii, line in enumerate(left_assignments):
    text += f'        {line}\n'
text += '        \n'

text += '    } else if (k > k2-1) {\n'
for ii, line in enumerate(right_assignments):
    text += f'        {line}\n'
text += '        \n'

text += '    } else {\n'
for ii, line in enumerate(interior_assignments):
    text += f'        {line}\n'
text += '    }\n'
text += '}\n'

f = open(f'{dirname}diffeq.cpp', 'w')
f.write(text)
f.close()

#############################
### WRITE MAIN.CPP FILE ###
#############################
text = ''
text += '#include <cmath>\n'
text += '#include <vector>\n'
text += '#include <iostream>\n'
text += '#include <fstream>\n'
text += '#include <chrono>\n'
text += '\n'
text += '#include "nr3.h"\n'
text += '#include "solvde.h"\n'
text += '\n'
text += '#include "diffeq.hpp"\n'
text += '\n'
text += 'using namespace std;\n'
text += '\n'

text += 'Int main(Int argc, const char * argv[]) {\n'
text += f'    Int N_grid = {N_grid};\n'
text += '    \n'

# Add for-loops
indent = '    '

for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Int':
        text += f'{indent}for ({paramdict_type[param]} {param.name}={paramdict_min[param]}; {param.name}<{paramdict_max[param]}+1; {param.name}++)' + ' {\n'
        indent += '    '

text += f'{indent}auto full_start = chrono::high_resolution_clock::now();\n'
text += f'{indent}float full_elapsed;\n'
text += f'{indent}\n'
text += f'{indent}VecDoub x_grid(N_grid);\n'
text += f'{indent}\n'

text += f'{indent}Doub h = ({x_max} - {x_min}) / ((Doub)N_grid-1.);\n'
text += f'{indent}for (Int ii=0; ii<N_grid; ii++)' + ' {\n'
text += f'{indent}    x_grid[ii] = (Doub)ii * h + {x_min};\n'
text += f'{indent}' + '}\n'

text += f'{indent}\n'

text += f'{indent}Int run_itmax = {run_itmax};\n'
text += f'{indent}Doub run_conv = {run_conv};\n'
text += f'{indent}Doub run_slowc = {run_slowc};\n'
text += f'{indent}VecDoub run_scalv({len(run_scalv)});\n'

for ii, val in enumerate(run_scalv):
    text += f'{indent}run_scalv[{ii}] = {val};\n'

text += f'{indent}VecInt run_indexv({len(Unknowns)});\n'

for ii in range(len(Unknowns)):
    text += f'{indent}run_indexv[{ii}] = {ii};\n'

text += f'{indent}Int run_NB = {len(Boundary_left)};\n'
text += f'{indent}\n'

text += f'{indent}MatDoub y_vec({len(Unknowns)},N_grid);'
text += f'{indent}\n'

for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Doub':
        text += f'{indent}Int {param.name}_points = {paramdict_points[param]};\n'
        text += f'{indent}VecDoub {param.name}_grid({param.name}_points);\n'
        text += f'{indent}for (Int ii=0; ii<{param.name}_points; ii++) ' + '{\n'
        text += f'{indent}    {param.name}_grid[ii] = ({paramdict_max[param]} - {paramdict_min[param]}) * ((Doub)ii / ((Doub){param.name}_points - 1.)) + {paramdict_min[param]};\n'
        text += f'{indent}' + '}\n'
        text += f'{indent}\n'

# File-writing
outname = f'"{outdir}'
for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Int':
        outname += f'{param.name}" + to_string({param.name}) + "_'

outname = outname[:-1] + '.txt"'

text += f'{indent}ofstream outfile({outname});\n'
text += f'{indent}outfile.precision({precision});\n'

first_cols = ' '.join([param.name for ii, param in enumerate(paramlist)])
first_cols += ' conv elapsed'

for ii, Unknown in enumerate(Unknowns):
    if Unknowns_flat[Unknown]:
        first_cols += f' {Unknown.name}'

text += f'{indent}outfile << "{first_cols}";\n'

for ii, Unknown in enumerate(Unknowns):
    if not Unknowns_flat[Unknown]:
        text += f'{indent}for (Int ii=0; ii<N_grid; ii++) ' + '{\n'
        text += f'{indent}    outfile << " {Unknown.name}" << to_string(ii);\n'
        text += f'{indent}' + '}\n'

text += f"{indent}outfile << '" + r'\n' + "';\n"
text += f'{indent}outfile.close();\n'
text += f'{indent}\n'

text += f'{indent}Int convtest;\n'
text += f'{indent}float elapsed;\n'
text += f'{indent}float thresh = {thresh};\n'
text += f'{indent}\n'

for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Doub':
        text += f'{indent}for (Int {param.name}_ii=0; {param.name}_ii<{param.name}_points; {param.name}_ii++)' + ' {\n'
        indent += '    '

text += f'{indent}auto start = chrono::high_resolution_clock::now();\n'
text += f'{indent}convtest = 1;\n'
text += f'{indent}\n'

text += f'{indent}Diffeq diffeq('
for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Int':
        text += f'{param}, '
    if paramdict_type[param] == 'Doub':
        text += f'{param}_grid[{param}_ii], '
    
text += 'x_grid);\n'
text += f'{indent}\n'

text += f'{indent}for (Int ii=0; ii<N_grid; ii++)' + ' {\n'
for jj, Unknown in enumerate(Unknowns):
    if bumps[jj] != 0:
        text += f'{indent}    y_vec[{jj}][ii] = {bumps[jj]};\n'
text += f'{indent}' + '}\n'
text += f'{indent}\n'

# INITIAL GUESSES (do in each loop)
text += f'{indent}for (Int ii=0; ii<N_grid; ii++)' + ' {\n'
for jj, guess in enumerate(Guesses):
    text += f'{indent}    y_vec[{jj}][ii] = {guess};\n'
text += f'{indent}' + '}\n'
text += f'{indent}\n'

text += f'{indent}Solvde<Diffeq> run_solvde(run_itmax, run_conv, run_slowc, run_scalv, run_indexv, run_NB, y_vec, diffeq, convtest);\n'
text += f'{indent}\n'

text += f'{indent}auto stop = chrono::high_resolution_clock::now();\n'
text += f'{indent}elapsed = (float)((stop - start).count());'
text += f'{indent}\n'

text += f'{indent}ofstream outfile({outname}, std::ios_base::app);\n'
text += f'{indent}outfile.precision({precision});\n'
text += f'{indent}\n'

text += f'{indent}outfile'
for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Int':
        text += f' << to_string({param.name}) << " "'
    if paramdict_type[param] == 'Doub':
        text += f' << to_string({param.name}_grid[{param.name}_ii]) << " "'

text += ' << convtest << " " << elapsed'

for ii, Unknown in enumerate(Unknowns):
    if Unknowns_flat[Unknown]:
        text += f' << " " << to_string(y_vec[{ii}][0])'
text += ';\n'

for ii, Unknown in enumerate(Unknowns):
    if not Unknowns_flat[Unknown]:
        text += f'{indent}for (Int jj=0; jj<N_grid; jj++) ' + '{\n'
        text += f'{indent}    outfile << " " << to_string(y_vec[{ii}][jj]);\n'
        text += indent + '}\n'

text += f"{indent}outfile << '" + r'\n' + "';\n"
text += f'{indent}outfile.close();\n'

for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Doub':
        indent = indent[4:]
        text += indent + '}\n'

text += f'{indent}\n'
text += f'{indent}auto full_stop = chrono::high_resolution_clock::now();\n'
text += f'{indent}full_elapsed = (float)((full_stop - full_start).count());\n'

for ii, param in enumerate(paramlist):
    if paramdict_type[param] == 'Int':
        indent = indent[4:]
        text += indent + '}\n'

text += '    \n'
text += '    return 0;\n'

text += '}\n'

f = open(f'{dirname}main.cpp', 'w')
f.write(text)
f.close()
