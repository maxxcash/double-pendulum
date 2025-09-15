import sympy as sp
import numpy as np
import symengine as se


class Dynamics_IDP:
    def __init__(self, I1c, I2c, I3c, m1, m2, m3, d1, d2, d3, l1, l2, l3, g):
        self.I1c, self.I2c, self.I3c = I1c, I2c, I3c
        self.m1, self.m2, self.m3 = m1, m2, m3
        self.d1, self.d2, self.d3 = d1, d2, d3
        self.l1, self.l2, self.l3 = l1, l2, l3
        self.g = g
        self.setup_symbols()
        self.I_func = None
        self.h_func = None
        self.V_pend_func = None
        self.gamma_func = None
        self.delta_func = None
        self._I = None
        self._h = None
        self._V_pend = None
        self._ext_f = None
        self._theta_d_d = None
        self._gamma = None
        self._delta = None

    def setup_symbols(self):
        self.t = sp.Symbol('t', real=True)
        self.theta1 = sp.Function('theta1')(self.t)
        self.theta2 = sp.Function('theta2')(self.t)
        self.theta3 = sp.Function('theta3')(self.t)
        self.tau = sp.Symbol('tau', real=True)
        self.x1, self.x2, self.x3, self.x4, self.x5, self.x6 = sp.symbols('x1 x2 x3 x4 x5 x6', real=True)




    def vec_cross(self, a):
        return sp.Matrix([
            [0, -a[2], a[1]],
            [a[2], 0, -a[0]],
            [-a[1], a[0], 0]
        ])

    def compute_I_and_h(self):
        if self._I is None or self._h is None or self._ext_f is None:
            g_f = sp.Matrix([0, 0, -self.g])
            Id = sp.eye(3)
            Id_1 = sp.eye(6)
            ZO = sp.zeros(3, 3)
            ZO_1 = sp.zeros(6, 6)
            ZO_c = sp.zeros(6, 1)
            ZO_c2 = sp.zeros(3, 1)

            I1 = self.I1c + ((self.m1 * self.d1 * self.d1) * Id) - (self.m1 * np.outer([self.d1, 0, 0], [self.d1, 0, 0]))
            I2 = self.I2c + ((self.m2 * self.d2 * self.d2) * Id) - (self.m2 * np.outer([self.d2, 0, 0], [self.d2, 0, 0]))
            I3 = self.I3c + ((self.m3 * self.d3 * self.d3) * Id) - (self.m3 * np.outer([self.d3, 0, 0], [self.d3, 0, 0]))

            ang = sp.Symbol('ang')
            Rz = sp.Matrix([
                [sp.cos(ang), -sp.sin(ang), 0],
                [sp.sin(ang), sp.cos(ang), 0],
                [0, 0, 1]
            ])
            Ry = sp.Matrix([
                [sp.cos(ang), 0, sp.sin(ang)],
                [0, 1, 0],
                [-sp.sin(ang), 0, sp.cos(ang)]
            ])

            T2_1 = Rz.subs(ang, self.theta1)
            T3_1 = T2_1 * Ry.subs(ang, sp.pi/2)
            T4_1 = T3_1 * Rz.subs(ang, self.theta2)
            T5_1 = T4_1 * Rz.subs(ang, self.theta3)

            d_1 = T2_1 * sp.Matrix([self.d1, 0, 0])
            d_2 = T4_1 * sp.Matrix([self.d2, 0, 0])
            d_3 = T5_1 * sp.Matrix([self.d3, 0, 0])

            e1 = sp.Matrix([0, 0, 1])
            e2 = sp.Matrix([0, 0, 1])
            e3 = sp.Matrix([0, 0, 1])
            p1 = sp.Matrix.vstack(e1, ZO_c2)
            p2 = sp.Matrix.vstack(T4_1 * e2, ZO_c2)
            p3 = sp.Matrix.vstack(T5_1 * e3, ZO_c2)
  
            a2_1 = T2_1 * sp.Matrix([-self.l1, 0, 0])
            a3_2 = T4_1 * sp.Matrix([-self.l2, 0, 0])
            a3_1 = a3_2 + a2_1
            A21 = sp.Matrix.vstack(sp.Matrix.hstack(Id, ZO), sp.Matrix.hstack(self.vec_cross(a2_1), Id))
            A31 = sp.Matrix.vstack(sp.Matrix.hstack(Id, ZO), sp.Matrix.hstack(self.vec_cross(a3_1), Id))
            A32 = sp.Matrix.vstack(sp.Matrix.hstack(Id, ZO), sp.Matrix.hstack(self.vec_cross(a3_2), Id))
            Nl = sp.Matrix.vstack(
                sp.Matrix.hstack(Id_1, ZO_1, ZO_1),
                sp.Matrix.hstack(A21, Id_1, ZO_1),
                sp.Matrix.hstack(A31, A32, Id_1)
            )
            Nd = sp.Matrix.vstack(
                sp.Matrix.hstack(p1, ZO_c, ZO_c),
                sp.Matrix.hstack(ZO_c, p2, ZO_c),
                sp.Matrix.hstack(ZO_c, ZO_c, p3)
            )

            I1_f = T2_1 * I1 * T2_1.T
            I2_f = T4_1 * I2 * T4_1.T
            I3_f = T5_1 * I3 * T5_1.T

            M1 = sp.Matrix.vstack(
                sp.Matrix.hstack(I1_f, self.m1 * self.vec_cross(d_1)),
                sp.Matrix.hstack(-self.m1 * self.vec_cross(d_1), self.m1 * Id)
            )
            M2 = sp.Matrix.vstack(
                sp.Matrix.hstack(I2_f, self.m2 * self.vec_cross(d_2)),
                sp.Matrix.hstack(-self.m2 * self.vec_cross(d_2), self.m2 * Id)
            )
            M3 = sp.Matrix.vstack(
                sp.Matrix.hstack(I3_f, self.m3 * self.vec_cross(d_3)),
                sp.Matrix.hstack(-self.m3 * self.vec_cross(d_3), self.m3 * Id)
            )
            M = sp.Matrix.vstack(
                sp.Matrix.hstack(M1, ZO_1, ZO_1),
                sp.Matrix.hstack(ZO_1, M2, ZO_1),
                sp.Matrix.hstack(ZO_1, ZO_1, M3)
            )

            N = Nl * Nd
            N_d = N.diff(self.t)

            ang1 = sp.Matrix([0, 0, self.theta1]).diff(self.t)
            ang2 = ang1 + (T4_1 * sp.Matrix([0, 0, self.theta2]).diff(self.t))
            ang3 = ang2 + (T5_1 * sp.Matrix([0, 0, self.theta3]).diff(self.t))

            w1 = sp.Matrix.vstack(
                sp.Matrix.hstack(self.vec_cross(ang1), ZO),
                sp.Matrix.hstack(ZO, self.vec_cross(ang1))
            )
            w2 = sp.Matrix.vstack(
                sp.Matrix.hstack(self.vec_cross(ang2), ZO),
                sp.Matrix.hstack(ZO, self.vec_cross(ang2))
            )
            w3 = sp.Matrix.vstack(
                sp.Matrix.hstack(self.vec_cross(ang3), ZO),
                sp.Matrix.hstack(ZO, self.vec_cross(ang3))
            )
            W = sp.Matrix.vstack(
                sp.Matrix.hstack(w1, ZO_1, ZO_1),
                sp.Matrix.hstack(ZO_1, w2, ZO_1),
                sp.Matrix.hstack(ZO_1, ZO_1, w3)
            )

            E1 = sp.Matrix.vstack(sp.Matrix.hstack(Id, ZO), sp.Matrix.hstack(ZO, ZO))
            E = sp.Matrix.vstack(
                sp.Matrix.hstack(E1, ZO_1, ZO_1),
                sp.Matrix.hstack(ZO_1, E1, ZO_1),
                sp.Matrix.hstack(ZO_1, ZO_1, E1)
            )

            joint_r = sp.Matrix([self.theta1.diff(self.t), self.theta2.diff(self.t), self.theta3.diff(self.t)])
            twist = N * joint_r

            W1_e = sp.Matrix.vstack(d_1.cross(self.m1 * g_f), self.m1 * g_f)
            W2_e = sp.Matrix.vstack(d_2.cross(self.m2 * g_f), self.m2 * g_f)
            W3_e = sp.Matrix.vstack(d_3.cross(self.m3 * g_f), self.m3 * g_f)
            Wrench = sp.Matrix.vstack(W1_e, W2_e, W3_e)
            ext_f = sp.Matrix([self.tau, 0, 0])

            # *** Convert to SymEngine Matrices for faster calculations ***
            N_se = se.Matrix(N)
            N_d_se = se.Matrix(N_d)
            M_se = se.Matrix(M)
            W_se = se.Matrix(W)
            E_se = se.Matrix(E)
            joint_r_se = se.Matrix(joint_r)
            Wrench_se = se.Matrix(Wrench)
            ext_f_se = se.Matrix(ext_f)

            # *** Perform Calculations using SymEngine ***
            I_se = N_se.T * M_se * N_se
            h_se = (N_se.T * M_se * N_d_se * joint_r_se) + (N_se.T * W_se * M_se * E_se * twist) - (N_se.T * Wrench_se)

            # # Define the substitution dictionary
            subs_dict = {
                self.theta1.diff(self.t): self.x4,
                self.theta2.diff(self.t): self.x5,
                self.theta3.diff(self.t): self.x6,
                self.theta1: self.x1, 
                self.theta2: self.x2,
                self.theta3: self.x3
            }

            # Apply the substitutions
            I_sub = I_se.replace(subs_dict)
            h_sub = h_se.replace(subs_dict)
            ext_f_sub = ext_f_se.replace(subs_dict)

            #Convert the final result to sympy matrix
            # I = sp.Matrix(I_sub)
            # h = sp.Matrix(h_sub)
            # ext_f = sp.Matrix(ext_f_sub)
            I = (I_sub)
            h = (h_sub)
            ext_f = (ext_f_sub)

            #I = sp.simplify(N.T * M * N)
            #h = sp.simplify((N.T * M * N_d * joint_r) + (N.T * W * M * E * twist) - (N.T * Wrench))
            # I = N.T * M * N
            # h = (N.T * M * N_d * joint_r) + (N.T * W * M * E * twist) - (N.T * Wrench)

            self._I = I
            self._h = h
            self._ext_f = ext_f

        return self._I, self._h, self._ext_f

    def compute_V_pend(self):
        if self._V_pend is None:

            ang = sp.Symbol('ang')
            Rz = sp.Matrix([
                [sp.cos(ang), -sp.sin(ang), 0],
                [sp.sin(ang), sp.cos(ang), 0],
                [0, 0, 1]
            ])
            Ry = sp.Matrix([
                [sp.cos(ang), 0, sp.sin(ang)],
                [0, 1, 0],
                [-sp.sin(ang), 0, sp.cos(ang)]
            ])

            T2_1 = Rz.subs(ang, self.theta1)
            T3_1 = T2_1 * Ry.subs(ang, sp.pi/2)
            T4_1 = T3_1 * Rz.subs(ang, self.theta2)
            T5_1 = T4_1 * Rz.subs(ang, self.theta3)

            d_1 = T2_1 * sp.Matrix([self.d1, 0, 0])
            d_2 = T4_1 * sp.Matrix([self.d2, 0, 0])
            d_3 = T5_1 * sp.Matrix([self.d3, 0, 0])

            ang1 = sp.Matrix([0, 0, self.theta1]).diff(self.t)
            ang2 = ang1 + (T4_1 * sp.Matrix([0, 0, self.theta2]).diff(self.t))
            ang3 = ang2 + (T5_1 * sp.Matrix([0, 0, self.theta3]).diff(self.t))

            g_f = sp.Matrix([0, 0, -self.g])

            l_1 = T2_1 * sp.Matrix([self.l1, 0, 0])
            v2_mag = ((-self.vec_cross(ang1) * l_1) + (-self.vec_cross(ang2) * d_2)).norm()
            KE2t = 0.5 * self.m2 * (v2_mag ** 2)
            KE2r = 0.5 * (ang2.T * (T4_1.T * self.I2c * T4_1) * ang2)[0]
            KE2 = KE2r + KE2t
            V2 = - (self.m2 * g_f).dot(d_2 + l_1)
            E2 = KE2 + V2

            subs_dict = {
                self.theta1.diff(self.t): self.x4,
                self.theta2.diff(self.t): self.x5,
                self.theta3.diff(self.t): self.x6,
                self.theta1: self.x1, 
                self.theta2: self.x2,
                self.theta3: self.x3
            }

            self._V_pend = E2.subs(subs_dict)

        return self._V_pend
    
    def compute_theta_d_d(self):
            if self._theta_d_d is None:
                I, h, ext_f = self.compute_I_and_h()
                #print(I[0:2,0:2])
                I_inv = I[0:2,0:2].inv()
                theta_d_d_inter = I_inv * (ext_f[0:2,0] - h[0:2,0])
                zero_row = se.Matrix([[0]])
                self._theta_d_d = theta_d_d_inter.row_insert(theta_d_d_inter.nrows(),zero_row)
            return self._theta_d_d

    def compute_gamma(self):
        if self._gamma is None:
            theta_d_d = self.compute_theta_d_d()
            self._gamma = se.diff(theta_d_d[0], self.tau)
        return self._gamma

    def compute_delta(self):
        if self._delta is None:
            theta_d_d = self.compute_theta_d_d()
            gamma = self.compute_gamma()
            self._delta = theta_d_d[0] - (self.tau * gamma)
        return self._delta

    def get_I_func(self):
        if self.I_func is None:
            I, _, _ = self.compute_I_and_h()
            I = sp.Matrix(I)
            # vars = (self.theta1, self.theta2, self.theta3, 
            #         self.theta1.diff(self.t), self.theta2.diff(self.t), self.theta3.diff(self.t))
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6)
            self.I_func = sp.lambdify(vars, I, modules=['numpy'])
        return self.I_func

    def get_h_func(self):
        if self.h_func is None:
            _, h, _ = self.compute_I_and_h()
            h = sp.Matrix(h)
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6)
            self.h_func = sp.lambdify(vars, h, modules=['numpy'])
        return self.h_func

    def get_gamma_func(self):
        if self.gamma_func is None:
            gamma = self.compute_gamma()   
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6,self.tau)
            self.gamma_func = sp.lambdify(vars, gamma, modules=['numpy'])
        return self.gamma_func

    def get_delta_func(self):
        if self.delta_func is None:
            delta = self.compute_delta()
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6,self.tau)
            self.delta_func = sp.lambdify(vars, delta, modules=['numpy'])
        return self.delta_func

    def get_V_pend_func(self):
        if self.V_pend_func is None:
            V_pend = self.compute_V_pend()
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6)
            self.V_pend_func = sp.lambdify(vars, V_pend, modules=['numpy'])
        return self.V_pend_func



# # Example values for link features (replace with your actual values)
# I1 = np.diag([21e-6, 133e-6, 141e-6])  # Inertia matrix for link 1
# I2 = np.diag([2e-6, 72e-6, 73e-6])     # Inertia matrix for link 2
# I3 = np.diag([0, 0, 0])     # Inertia matrix for link 3
# m1, m2, m3 = 0.069, 0.036, 0       # Masses of links
# d1, d2, d3 = 0.02, 0.042, 0        # Center of mass distances
# l1, l2, l3 = 0.085, 0.134, 0       # Link lengths
# g = 9.81                               # Gravity

# state=[0, 0, 0, 0, 0, 0]
# desired_state=[0, np.pi, 0, 0, 0, 0]
# control_obj = Dynamics_IDP(I1, I2, I3, m1, m2, m3, d1, d2, d3, l1, l2, l3, g)
    
# I_func = control_obj.get_I_func()
# print('Process-1 done')
# h_func = control_obj.get_h_func()
# print('Process-2 done')
# gamma_func = control_obj.get_gamma_func()
# print('Process-3 done')
# delta_func = control_obj.get_delta_func()
# print('Process-4 done')

