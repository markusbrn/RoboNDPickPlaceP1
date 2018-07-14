from sympy import *

q1, q2, q3, q4, q5, q6, qG = symbols('q1:8')
d1, d2, d3, d4, d5, d6, dG = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
# Create Modified DH parameters
s = {alpha0:     0,   a0:      0,   d1:  0.75,
     alpha1: -pi/2,   a1:   0.35,   d2:     0,   q2: q2-pi/2,
     alpha2:     0,   a2:   1.25,   d3:     0,
     alpha3: -pi/2,   a3: -0.054,   d4:  1.50,
     alpha4:  pi/2,   a4:      0,   d5:     0,
     alpha5: -pi/2,   a5:      0,   d6:     0,
     alpha6:     0,   a6:      0,   dG: 0.303,   qG: 0}
# Define Modified DH Transformation matrix
def Matrix_DH(alpha, a, d, q):
    return Matrix([[           cos(q),           -sin(q),           0,             a],
                   [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                   [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                   [                0,                 0,           0,             1]])
# Create individual transformation matrices
T_0_1 = simplify(Matrix_DH(alpha0, a0, d1, q1).subs(s))
#print(T_0_1)
T_1_2 = simplify(Matrix_DH(alpha1, a1, d2, q2).subs(s))
#print(T_1_2)
T_2_3 = simplify(Matrix_DH(alpha2, a2, d3, q3).subs(s))
#print(T_2_3)
T_3_4 = simplify(Matrix_DH(alpha3, a3, d4, q4).subs(s))
#print(T_3_4)
T_4_5 = simplify(Matrix_DH(alpha4, a4, d5, q5).subs(s))
#print(T_4_5)
T_5_6 = simplify(Matrix_DH(alpha5, a5, d6, q6).subs(s))
#print(T_5_6)
T_6_G = simplify(Matrix_DH(alpha6, a6, dG, qG).subs(s))
#print(T_6_G)
T_0_G = T_0_1*T_1_2*T_2_3*T_3_4*T_4_5*T_5_6*T_6_G
#print(pretty(simplify(T_0_G)))

R_3_G = pretty(simplify(T_3_4[0:3,0:3] * T_4_5[0:3,0:3] * T_5_6[0:3,0:3] * T_6_G[0:3,0:3]))
#print(R_3_G)

R_4_G = pretty(simplify(T_4_5[0:3,0:3] * T_5_6[0:3,0:3] * T_6_G[0:3,0:3]))
print(R_4_G)

R_5_G = pretty(simplify(T_5_6[0:3,0:3] * T_6_G[0:3,0:3]))
print(R_5_G)
