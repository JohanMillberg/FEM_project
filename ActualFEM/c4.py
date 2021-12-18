from dolfin import *
import random
import matplotlib.pyplot as plt

mesh = Mesh("sphere2.xml")

V = VectorFunctionSpace(mesh, 'P', 1, dim=2)

T = 20
dt = 0.1
num_steps = int(T/dt)
rho = 2
alpha_1 = 0.01
alpha_2 = 0.005
cf = 0.024
ck = 0.055
R = 0.5
r = 0.3


class InitialConditions(UserExpression):

    def eval(self, values, x):
        if (R-(x[0]**2+x[1]**2)**0.5)**2 + x[2]**2 <= r**2:
            values[0] = rho*1/2*random.uniform(0,1)
            values[1] = rho*(1+1/2*random.uniform(0,1))

        else:
            values[0] = 0.0
            values[1] = 0.0

    def value_shape(self):
        return (2,)

indata = InitialConditions(degree=2)
u0 = Function(V)
u0 = interpolate(indata, V)
u_n, v_n = split(u0)

U = TrialFunction(V)
u, v = split(U)

h, g = TestFunctions(V)

#Creates constants
dt = Constant(dt)
alpha_1 = Constant(alpha_1)
alpha_2 = Constant(alpha_2)
cf = Constant(cf)
ck = Constant(ck)

#Creates linear and bilinear form

#F1 = ((u-u_n)/dt*h + alpha_1*inner(grad(u),grad(h))+cf*u*h - h*(-u_n*v_n*v_n+cf))*dx
#Using Crank-Nicolson:
F1 = ((u-u_n)/dt*h + alpha_1/2*inner(grad(u+u_n),grad(h))+cf/2*(u+u_n)*h - h*(-u_n*v_n*v_n+cf))*dx

#F2 = ((v-v_n)/dt*g + alpha_2*inner(grad(v),grad(g))+(cf+ck)*v*g-g*(u_n*v_n*v_n))*dx
#Using Crank-Nicolson:
F2 = ((v-v_n)/dt*g + alpha_2/2*inner(grad(v+v_n),grad(g))+(cf+ck)/2*(v+v_n)*g-g*(u_n*v_n*v_n))*dx

F = F1 + F2

a, L = lhs(F), rhs(F)

U = Function(V)

vtkfile_u = File("solutions_u3D/solution_u.pvd")
vtkfile_v = File("solutions_v3D/solution_v.pvd")

u_initial = interpolate(indata, V)

mass_u = []
mass_v = []
time_vec = []

t = 0
for n in range(num_steps+1):

    t += dt
    time_vec.append(t)
    solve(a == L, U)
    u_, v_ = U.split()
    vtkfile_u << (u_, t)
    vtkfile_v << (v_, t)
    u0.assign(U)

    M_u = (u_initial.split()[0] - u_)*dx
    M_v = (u_initial.split()[1] - v_)*dx

    #M_u = (U.split()[0] - u_initial.split()[0])*dx
    #M_v = (U.split()[1] - u_initial.split()[1])*dx

    mass_u.append(assemble(M_u))
    mass_v.append(assemble(M_v))

fig, ax = plt.subplots()
plt.plot(time_vec,mass_u,label='u')
plt.plot(time_vec,mass_v,label='v')
ax.set_xlabel('Time')
ax.set_ylabel('Mass Loss')
ax.set_title('Mass loss over time of u and v')
plt.legend()
plt.show()