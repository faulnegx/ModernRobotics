import numpy as np
# Skip to Line 161. Everything before is all copied from modern robot library.
def NearZero(z):
    return abs(z) < 1e-6

def Normalize(V):
    return V / np.linalg.norm(V)

def so3ToVec(so3mat):
    return np.array([so3mat[2][1], so3mat[0][2], so3mat[1][0]])

def AxisAng3(expc3):
    return (Normalize(expc3), np.linalg.norm(expc3))

def MatrixExp3(so3mat):
    omgtheta = so3ToVec(so3mat)
    if NearZero(np.linalg.norm(omgtheta)):
        return np.eye(3)
    else:
        theta = AxisAng3(omgtheta)[1]
        omgmat = so3mat / theta
        return np.eye(3) + np.sin(theta) * omgmat \
               + (1 - np.cos(theta)) * np.dot(omgmat, omgmat)

def VecToso3(omg):
    return np.array([[0,      -omg[2],  omg[1]],
                     [omg[2],       0, -omg[0]],
                     [-omg[1], omg[0],       0]])

def VecTose3(V):
    return np.r_[np.c_[VecToso3([V[0], V[1], V[2]]), [V[3], V[4], V[5]]],
                 np.zeros((1, 4))]

def MatrixExp6(se3mat):
    se3mat = np.array(se3mat)
    omgtheta = so3ToVec(se3mat[0: 3, 0: 3])
    if NearZero(np.linalg.norm(omgtheta)):
        return np.r_[np.c_[np.eye(3), se3mat[0: 3, 3]], [[0, 0, 0, 1]]]
    else:
        theta = AxisAng3(omgtheta)[1]
        omgmat = se3mat[0: 3, 0: 3] / theta
        return np.r_[np.c_[MatrixExp3(se3mat[0: 3, 0: 3]),
                           np.dot(np.eye(3) * theta \
                                  + (1 - np.cos(theta)) * omgmat \
                                  + (theta - np.sin(theta)) \
                                    * np.dot(omgmat,omgmat),
                                  se3mat[0: 3, 3]) / theta],
                     [[0, 0, 0, 1]]]

def FKinSpace(M, Slist, thetalist):
    T = np.array(M)
    for i in range(len(thetalist) - 1, -1, -1):
        T = np.dot(MatrixExp6(VecTose3(np.array(Slist)[:, i] \
                                       * thetalist[i])), T)
    return T


def TransToRp(T):
    T = np.array(T)
    return T[0: 3, 0: 3], T[0: 3, 3]


def Adjoint(T):
    R, p = TransToRp(T)
    return np.r_[np.c_[R, np.zeros((3, 3))],
                 np.c_[np.dot(VecToso3(p), R), R]]

def se3ToVec(se3mat):
    return np.r_[[se3mat[2][1], se3mat[0][2], se3mat[1][0]],
                 [se3mat[0][3], se3mat[1][3], se3mat[2][3]]]

def MatrixLog3(R):
    acosinput = (np.trace(R) - 1) / 2.0
    if acosinput >= 1:
        return np.zeros((3, 3))
    elif acosinput <= -1:
        if not NearZero(1 + R[2][2]):
            omg = (1.0 / np.sqrt(2 * (1 + R[2][2]))) \
                  * np.array([R[0][2], R[1][2], 1 + R[2][2]])
        elif not NearZero(1 + R[1][1]):
            omg = (1.0 / np.sqrt(2 * (1 + R[1][1]))) \
                  * np.array([R[0][1], 1 + R[1][1], R[2][1]])
        else:
            omg = (1.0 / np.sqrt(2 * (1 + R[0][0]))) \
                  * np.array([1 + R[0][0], R[1][0], R[2][0]])
        return VecToso3(np.pi * omg)
    else:
        theta = np.arccos(acosinput)
        return theta / 2.0 / np.sin(theta) * (R - np.array(R).T)


def MatrixLog6(T):
    R, p = TransToRp(T)
    omgmat = MatrixLog3(R)
    if np.array_equal(omgmat, np.zeros((3, 3))):
        return np.r_[np.c_[np.zeros((3, 3)),
                           [T[0][3], T[1][3], T[2][3]]],
                     [[0, 0, 0, 0]]]
    else:
        theta = np.arccos((np.trace(R) - 1) / 2.0)
        return np.r_[np.c_[omgmat,
                           np.dot(np.eye(3) - omgmat / 2.0 \
                           + (1.0 / theta - 1.0 / np.tan(theta / 2.0) / 2) \
                              * np.dot(omgmat,omgmat) / theta,[T[0][3],
                                                               T[1][3],
                                                               T[2][3]])],
                     [[0, 0, 0, 0]]]


def JacobianSpace(Slist, thetalist):
    Js = np.array(Slist).copy().astype(np.float)
    T = np.eye(4)
    for i in range(1, len(thetalist)):
        T = np.dot(T, MatrixExp6(VecTose3(np.array(Slist)[:, i - 1] \
                                * thetalist[i - 1])))
        Js[:, i] = np.dot(Adjoint(T), np.array(Slist)[:, i])
    return Js

def JacobianBody(Blist, thetalist):
    Jb = np.array(Blist).copy().astype(np.float)
    T = np.eye(4)
    for i in range(len(thetalist) - 2, -1, -1):
        T = np.dot(T,MatrixExp6(VecTose3(np.array(Blist)[:, i + 1] \
                                         * -thetalist[i + 1])))
        Jb[:, i] = np.dot(Adjoint(T), np.array(Blist)[:, i])
    return Jb


def TransInv(T):
    R, p = TransToRp(T)
    Rt = np.array(R).T
    return np.r_[np.c_[Rt, -np.dot(Rt, p)], [[0, 0, 0, 1]]]

def FKinBody(M, Blist, thetalist):
    T = np.array(M)
    for i in range(len(thetalist)):
        T = np.dot(T, MatrixExp6(VecTose3(np.array(Blist)[:, i] \
                                          * thetalist[i])))
    return T

def IKinBody(Blist, M, T, thetalist0, eomg, ev):
    thetalist = np.array(thetalist0).copy()
    i = 0
    maxiterations = 20
    Vb = se3ToVec(MatrixLog6(np.dot(TransInv(FKinBody(M, Blist, \
                                                      thetalist)), T)))
    err = np.linalg.norm([Vb[0], Vb[1], Vb[2]]) > eomg \
          or np.linalg.norm([Vb[3], Vb[4], Vb[5]]) > ev
    while err and i < maxiterations:
        thetalist = thetalist \
                    + np.dot(np.linalg.pinv(JacobianBody(Blist, \
                                                         thetalist)), Vb)
        i = i + 1
        Vb \
        = se3ToVec(MatrixLog6(np.dot(TransInv(FKinBody(M, Blist, \
                                                       thetalist)), T)))
        err = np.linalg.norm([Vb[0], Vb[1], Vb[2]]) > eomg \
              or np.linalg.norm([Vb[3], Vb[4], Vb[5]]) > ev
    return (thetalist, not err)

# This starts custom code. 

def IKinBodyIterates(Blist, M, T, thetalist0, eomg, ev):
    """
    :param Blist: The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
    :param M: The home configuration of the end-effector
    :param T: The desired end-effector configuration Tsd
    :param thetalist0: An initial guess of joint angles that are close to
                       satisfying Tsd
    :param eomg: A small positive tolerance on the end-effector orientation
                 error. The returned joint angles must give an end-effector
                 orientation error less than eomg
    :param ev: A small positive tolerance on the end-effector linear position
               error. The returned joint angles must give an end-effector
               position error less than ev
    :return thetalist: Joint angles that achieve T within the specified
                       tolerances,
    :return success: A logical value where TRUE means that the function found
                     a solution and FALSE means that it ran through the set
                     number of maximum iterations without finding a solution
                     within the tolerances eomg and ev.
    Uses an iterative Newton-Raphson root-finding method.
    The maximum number of iterations before the algorithm is terminated has
    been hardcoded in as a variable called maxiterations. It is set to 20 at
    the start of the function, but can be changed if needed.

    Example based on Chapter 6.2.2"
    Inputs:
        Blist = np.transpose(np.array(
            [[0,0,1,0,2,0],
            [0,0,1,0,1,0]]
        ))
        T = np.array(
            [[-0.5,-0.866,0,0.366],
            [0.866,-0.5,0,1.366],
            [0,0,1,0],
            [0,0,0,1]]
        )
        M = np.array(
            [[1,0,0,2],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1]]
        )
        thetalist0 = np.array([0, np.pi/6])
    Call:
        IKinBodyIterates(Blist, M, T, thetalist0, 0.001, 0.001)
    returns:
        (array([0.52358947, 1.57082969]), True)
    """
    theta_n = np.array(thetalist0).copy()

    i = 0
    maxiterations = 20
    # Calculate current configuration
    Tsb = FKinBody(M, Blist, theta_n)
    # 
    Vb = se3ToVec(MatrixLog6(np.dot(TransInv(Tsb), T)))
    # Calculate angular and linear error magnitudes
    omega_b = np.linalg.norm([Vb[0], Vb[1], Vb[2]])
    lin_err = np.linalg.norm([Vb[3], Vb[4], Vb[5]])
    # report if error is small enough. If it is small enough, then err will be false.
    err = omega_b > eomg or lin_err > ev
    
    # print statements
    print "Iteration 0: "
    print "joint vector: {}".format(theta_n)
    print "SE(3) end-effector config: \n{}".format(Tsb)
    print "error twist V_b: \n{}".format(Vb)
    print "angular error magnitude ||omega_b||: {:0.3f}".format(omega_b)
    print "linear error magnitude ||v_b||: {:0.3f}".format(lin_err)
    print "\n"
    # Add the new theta into thetalist to store
    thetalist = np.array(theta_n).copy()

    while err and i < maxiterations:
        # Uses an iterative Newton-Raphson to find the next theta set to try.
        theta_n = theta_n + np.dot(np.linalg.pinv(JacobianBody(Blist, theta_n)),Vb)
        # Add the new theta into thetalist to store
        thetalist = np.vstack((thetalist, theta_n))
        i = i + 1
        Tsb = FKinBody(M, Blist, thetalist[-1])
        Vb = se3ToVec(MatrixLog6(np.dot(TransInv(Tsb), T)))
        # Calculate angular and linear error magnitudes
        omega_b = np.linalg.norm([Vb[0], Vb[1], Vb[2]])
        lin_err = np.linalg.norm([Vb[3], Vb[4], Vb[5]])
        
        err = omega_b > eomg or lin_err > ev

        print "Iteration {}: ".format(i)
        print "joint vector: {}".format(thetalist[-1])
        print "SE(3) end-effector config: \n{}".format(Tsb)
        print "error twist V_b: \n{}".format(Vb)
        print "angular error magnitude ||omega_b||: {:0.3f}".format(omega_b)
        print "linear error magnitude ||v_b||: {:0.3f}".format(lin_err)
        print "\n"
    # End of while loop
    # Print out thetalist and export to csv
    print "thetalist: \n{}".format(thetalist)
    np.savetxt("iterates.csv", thetalist, delimiter=",")
    print "iterates.csv saved to current directory. "
    # return the final theta value and report if we converge.
    return (theta_n, not err)

if __name__ == '__main__':
    np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
    # Example 4.5 of Chapter 4.1.2 (Figure 4.6)
    L1 = 425/1000.0
    L2 = 392/1000.0
    W1 = 109/1000.0
    W2 = 82/1000.0
    H1 = 89/1000.0
    H2 = 95/1000.0

    Blist = np.transpose(np.array([
        [0, 1,  0, W1+W2, 0,          L1+L2],
        [0, 0,  1, H2,    -1*(L1+L2), 0],
        [0, 0,  1, H2,    -1*L2,      0],
        [0, 0,  1, H2,    0,          0],
        [0, -1, 0, -1*W2, 0,          0],
        [0, 0,  1, 0,     0,          0]]
    ))
    T = np.array([
        [0,  1, 0,  -0.5],
        [0,  0, -1, 0.1],
        [-1, 0, 0,  0.1],
        [0,  0, 0,  1]]
    )

    M = np.array([
        [-1, 0, 0, L1+L2],
        [0,  0, 1, W1+W2],
        [0,  1, 0, H1-H2],
        [0,  0, 0, 1]]
    )
    thetalist0 = np.array([3, 0, np.pi/3, np.pi, np.pi/6, np.pi/4])
    print "Initial theta is {}.".format(thetalist0)
    print "Call 'print IKinBodyIterates(Blist, M, T, thetalist0, 0.001, 0.0001)' in script. \n"
    print IKinBodyIterates(Blist, M, T, thetalist0, 0.001, 0.0001)
