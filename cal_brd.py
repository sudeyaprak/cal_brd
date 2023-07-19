def cal_brd(eph, brd):
    """
    Computes the satellite coordinates in the Topocentric Radiant System (TRS) frame for a given ephemeris and broadcast
    ephemeris parameters using the given formulae. 

    Args:
        eph (float): The GPS time of interest (in seconds of day).
        brd (ndarray): A 7x4 array of broadcast ephemeris parameters for the satellite of interest. The rows represent
        the following parameters: (1) Crs, (2) delta_n, (3) M0, (4) Cuc, (5) e, (6) Cus, (7) sqrt(a), (8) toe, (9) Cic,
        (10) Omega0, (11) Cis, (12) i0, (13) Crc, (14) w, (15) Omega_dot, (16) i_dot.

    Returns:
        ndarray: A 3x1 array of the satellite's coordinates in the TRS frame at the given time.
    """
    
    # Extract the broadcast ephemeris parameters
    crs = brd[0,1] # Orbital Radius Correction (m)
    deln = brd[0,2] # Mean Motion Difference (rad/sec)
    mo = brd[0,3] # Mean Anomaly at Reference Epoch (rad)
    cuc = brd[1,0] # Latitude Argument Correction (rad)
    e = brd[1,1] # Eccentricity 
    cus = brd[1,2] # Latitude Argument Correction (rad)
    sqrt_a = brd[1,3] # Square Root of Semi-Major Axis
    toe = brd[2,0] # Ephemerides Reference Epoch in Seconds of Gps Week
    cic = brd[2,1] # Inclination Correction (rad)
    omega0 = brd[2,2] # Longitude of Ascending Note at the Beginning of the Week (rad)
    cis = brd[2,3] # Inclination Correction (rad)
    i0 = brd[3,0] # Inclination of Reference Epoch (rad)
    crc = brd[3,1] # Orbital Radius Correction (m)
    w = brd[3,2] # Argument of Perigee (rad)
    omegadot = brd[3,3] # Rate of Node's Right Assencion (rad/sec)
    idot = brd[4,0] # Rate of Inclination Angle (rad/sec)
    
    # Compute the time tk from the ephemerides reference epoch toe
    toe -= 3*(24*60*60) # seconds of gps week to seconds of day for 1 march wednesday
    t = eph
    tk = t - toe
    
    # Compute the mean anomaly for tk
    u = 3.986005e14 # Earth's gravitational constant (WGS84)
    mk = mo + ((np.sqrt(u)/(sqrt_a**3)) + deln)*tk
    
    # Solve (iteratively) the Kepler equation for the eccentric anomaly E_old
    E_new = mk
    E_old = 0
    while np.abs(E_new - E_old) > 1e-12:
        E_old = E_new
        E_new = E_old + (mk - E_old + e * math.sin(E_old)) / (1 - e * math.cos(E_old))
    
    # Compute the true anomaly vk
    vk = math.atan2(np.sqrt(1 - e**2) * math.sin(E_new), math.cos(E_new) - e)
    
    # Compute the argument of latitude uk
    uk = w + vk + cuc*math.cos(2*(w+vk)) + cus*math.sin(2*(w+vk))
    
    # Compute the radial distance rk
    a = sqrt_a**2
    rk = a*(1 - e*math.cos(E_new)) + crc*math.cos(2*(w+vk)) + crs*math.sin(2*(w+vk))
    
    # Compute the inclination ik
    ik = i0 + idot*tk + cic*math.cos(2*(w+vk)) + cis*math.sin(2*(w+vk))
    
    # Compute the longitude of ascending note
    we = 7.2921151467e-5 # Earth's rotation rate (WGS84) 
    lambdak = omega0 + (omegadot - we)*tk - (we*brd[2,0])
    
    # Define elementary rotation matrices and reflection matrices
    def getR1(alpha):
        return np.matrix([[1,0,0],
                          [0,np.cos(alpha),np.sin(alpha)],
                          [0,-np.sin(alpha),np.cos(alpha)]])
    
    def getR3(alpha):
        return np.matrix([[np.cos(alpha),np.sin(alpha),0],
                          [-np.sin(alpha),np.cos(alpha),0],
                          [0,0,1]])

    # Compute the coordinates in the TRS frame
    coor = getR3(-lambdak)*getR1(-ik)*getR3(-uk)*([[rk],[0], [0]])
    
    return coor
