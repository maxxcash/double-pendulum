import numpy as np
import sympy as sp
import symengine as se

from Dynamics_singlePend import Dynamics_IDP

class Controller(Dynamics_IDP):
    def __init__(self, I1, I2, I3, m1, m2, m3, d1, d2, d3, l1, l2, l3, g, state, desired_state):
        super().__init__(I1, I2, I3, m1, m2, m3, d1, d2, d3, l1, l2, l3, g)
        self.state = np.array(state)
        self.desired_state = np.array(desired_state)
        
        
        # Get functions for I, h, gamma, and delta
        self.I_func = self.get_I_func()
        self.h_func = self.get_h_func()
        self.gamma_func = self.get_gamma_func()
        self.delta_func = self.get_delta_func()
        self.V_pend_func = self.get_V_pend_func()

        # A and B for desired state
        self.f_x_x_func = None
        self.f_x_u_func = None
        self.A = None
        self.B = None
        self.K_gain = None
        self.t_switch = 0
        theta1_Zoff = 0
        theta2_Zoff = 0
        theta3_Zoff = 0
        self.zero_error = np.array([theta1_Zoff,theta2_Zoff,theta3_Zoff,0,0,0])

    def compute_torque(self, t, cur_state):

        self.state = cur_state
        theta1, theta2, theta3, dtheta1, dtheta2, dtheta3 = self.state
        theta1, theta2, theta3 = self.angle_constrain(theta1), self.angle_constrain(theta2), self.angle_constrain(theta3)
        theta1_des, theta2_des, theta3_des, dtheta1_des, dtheta2_des, dtheta3_des = self.desired_state

        # Compute numeric values for I, h, gamma, and delta
        I = self.I_func(theta1, theta2, theta3, dtheta1, dtheta2, dtheta3)
        h = self.h_func(theta1, theta2, theta3, dtheta1, dtheta2, dtheta3)
        gamma = self.gamma_func(theta1, theta2, theta3, dtheta1, dtheta2, dtheta3, 0)  # Assuming tau=0 for gamma
        delta = self.delta_func(theta1, theta2, theta3, dtheta1, dtheta2, dtheta3, 0)  # Assuming tau=0 for delta

        # Energy of pendulum
        V_p = self.V_pend_func(theta1, theta2, theta3, dtheta1, dtheta2, dtheta3)
        V_desired = 0.155*9.81*0.041


        # Define K_gain vector
        K_gain_rise = np.array([0.3162,0.8558])

        # Compute v
        if self.A is None or self.B is None:
            self.compute_A_B()
            # Removing x3(theta3) and x6(theta3_dot) for single inverted pendulum
            self.A = np.array([[self.A[0,0], self.A[0,1], self.A[0,3], self.A[0,4]],[self.A[1,0], self.A[1,1], self.A[1,3], self.A[1,4]],[self.A[3,0], self.A[3,1], self.A[3,3], self.A[3,4]],[self.A[4,0], self.A[4,1], self.A[4,3], self.A[4,4]]])
            self.B = np.array([[self.B[0,0]],[self.B[1,0]],[self.B[3,0]],[self.B[4,0]]])
        
        poles = [-2,-4,-2,-4]
        if self.K_gain is None:
            self.K_gain = self.ackermann(self.A,self.B,poles)
            

        states = np.array([theta1, theta2, theta3, dtheta1, dtheta2, dtheta3])
        desired_states = np.array([np.sign(theta1)*theta1_des, np.sign(theta2)*theta2_des, theta3_des, dtheta1_des, dtheta2_des, dtheta3_des])

        #Switching-based control
        if abs(states[1] - desired_states[1]) < 0.2 and abs(V_p - V_desired) < 0.061:
            #Stabilizing control
            print("Catch")
            mod_state = np.array([states[0],states[1],states[3],states[4]])
            mod_ds = np.array([desired_states[0],desired_states[1],desired_states[3],desired_states[4]])
            mod_K = np.array([self.K_gain[0],self.K_gain[1],self.K_gain[2],self.K_gain[3]])
            v = -np.dot(mod_K, (mod_state - mod_ds))
            #v = -np.dot(self.K_gain, (states - desired_states))  # Compute v as a scalar
            self.t_switch = t
        else:
            #tracking control - 1
            # print("0")
            # c = 1*(t%2)
            # a1, w1 = 1,c
            # const = (4*a1*w1*w1)
            # R = np.array([a1 * np.sin(w1 * t), a1 * w1 * np.cos(w1 * t)]) + np.array([self.zero_error[0],0])
            # err = np.array([states[0],states[3]]) - R
            # #print(err)
            # w = -np.dot(K_gain_rise,err)
            # #print('w = ',w)
            # v_d = w - (const*np.sin(w1*t)) + (2*a1*c*np.cos(w1*t))     # If w1 is a function of time
            # #v_d = w - (a1*w1*w1*np.sin(w1*t))       # For constant w1
            # v = (1 / gamma) * (v_d - delta)
            # v = v_d
            # # #print("gamma = ",gamma)
            # # #print("delta = ",delta)

            #tracking control - 2 (Based on paper by Yoshida) #This works, do not touch it
            c_o, eta = 0.9,1.2
            c = np.sqrt(9.81/0.042) #sqrt(g/d2)
            a_o, b_o = 1, 0.1
            print([abs(V_p - V_desired),(states[1] - desired_states[1])])

            if abs(V_p - V_desired) >= b_o:
                a1 = a_o*np.sign(V_p - V_desired)
            else:
                a1 = a_o*(V_p - V_desired)/b_o

            w1 = c
            g_wn = 0.417
            phi_wn = 1.571
            phi_t = w1*(t)
            const = a1/g_wn
            f1= (c/c_o)**2
            f2 = 2*np.sqrt(f1)*eta

            R = np.array([const * np.sin(phi_t - np.pi + phi_wn), const*w1*np.cos(phi_t - np.pi + phi_wn)]) + np.array([self.zero_error[0],0])
            v_d = (f1*(R[0] - states[0])) - (f2*(states[3]))
            v = v = (1 / gamma) * (v_d - delta)


        tau = v

        tau = max(min(tau, 0.40), -0.40)  # Clamp tau to [-0.40, 0.40]

        return tau

    def compute_A_B(self):
        if self.f_x_x_func is None or self.f_x_u_func is None:

            theta_dd = self._theta_d_d
            vars = (self.x1, self.x2, self.x3,self.x4,self.x5,self.x6,self.tau)
            f_x = [self.x4,self.x5,self.x6,theta_dd[0],theta_dd[1],theta_dd[2]]
            f_x = [self.x4,self.x5,self.x6,theta_dd[0],theta_dd[1],theta_dd[2]]
            f_x_x = se.Matrix([self.vec_diff(f_x,self.x1),self.vec_diff(f_x,self.x2),self.vec_diff(f_x,self.x3),self.vec_diff(f_x,self.x4),self.vec_diff(f_x,self.x5),self.vec_diff(f_x,self.x6)]).T
            f_x_u = se.Matrix(self.vec_diff(f_x,self.tau))

            theta1_des, theta2_des, theta3_des, dtheta1_des, dtheta2_des, dtheta3_des = self.desired_state
            subs_dict = {
                self.x1: theta1_des,
                self.x2: theta2_des,
                self.x3: theta3_des,
                self.x4: dtheta1_des,
                self.x5: dtheta2_des,
                self.x6: dtheta3_des,
                self.tau: 0 
            }
            self.A = np.array(f_x_x.subs(subs_dict))
            self.B = np.array(f_x_u.subs(subs_dict))
            self.f_x_x_func = f_x_x
            self.f_x_u_func = f_x_u
            
        
        else:

            theta1_des, theta2_des, theta3_des, dtheta1_des, dtheta2_des, dtheta3_des = self.desired_state
            subs_dict = {
                self.x1: theta1_des,
                self.x2: theta2_des,
                self.x3: theta3_des,
                self.x4: dtheta1_des,
                self.x5: dtheta2_des,
                self.x6: dtheta3_des,
                self.tau: 0 
            }
            self.A = np.array(f_x_x.subs(subs_dict))
            self.B = np.array(f_x_u.subs(subs_dict))

        return None

    def vec_diff(self,f,x):

        f_d = []
        for i in range(len(f)):
            f_d.append(se.diff(f[i],x))
        return f_d

    def ackermann(self, Ac, Bc, poles):

        n = Ac.shape[0]
        m = Bc.shape[1]

        # Forming the characteristic equation
        del_A = np.eye(n)
        for pole in poles:
            del_A = np.dot(del_A, (Ac - pole * np.eye(n)))

        # Forming the controllability matrix
        con = Bc
        iter = np.dot(Ac, Bc)
        for _ in range(n - 1):
            con = np.hstack((con, iter))
            iter = np.dot(Ac, iter)

        # Compute the state feedback gain K
        const = np.zeros((1, n))
        const[0, n - 1] = 1
        con = con.astype(np.float64)
        K = np.dot(const, np.linalg.inv(con)).dot(del_A)
        K = K.astype(np.float64)

        return K.flatten()

    def angle_constrain(self, theta):

        theta_con = (theta + np.pi) % (2*np.pi) - np.pi
        #theta_con = theta%(2*np.pi)

        return theta_con


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
# control_obj = Controller(I1, I2, I3, m1, m2, m3, d1, d2, d3, l1, l2, l3, g, state, desired_state)
# control_obj.compute_torque(0, state)

