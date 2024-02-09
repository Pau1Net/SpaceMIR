import math
from collections import namedtuple
#from units import Units

# Гравитационная постоянная (универсальная) (m^2/(kg*s^2))
G = 6.674e-11

# Радиус орбиты Земли или 1 астрономическая единица (AU) (m)
AU = 149.6e9

"""
class sun:
    def __init__(self, m, r_body):
        # Масса Солнца (kg)
        self.m = m
        # Радиус Солнца (m)
        self.r_body = r_body

class moon:
    def __init__(self, m, r_body, r_orbit):
        # Масса Луны (kg)
        self.m = m
        # Радиус Луны (m)
        self.r_body = r_body
        # Большая полуось (m)
        self.r_orbit = r_orbit

class planet:
    def __init__(self, m, r_body, r_orbit, b, phi_solar):
        # Масса планеты (выбор) (kg)
        self.m = m
        # Радиус планеты (m)
        self.r_body = r_body
        # Большая полуось (планета-Солнце) (m)
        self.r_orbit = r_orbit
        # Магнитное поле (T)
        self.b = b
        # Постоянный солнечный поток при ~1AU (W/m^2)
        self.phi_solar = phi_solar

# Солнце
sun = sun(m=1.989e30, r_body=695700e3)

# Земля
earth = planet(m=5.972e24, r_body=6371e3, r_orbit=1*AU, b=3e-5, phi_solar=1.3608e3)
# Луна
moon = moon(m=7.34e22, r_body=1.737e6, r_orbit=384748e3)

# Плутон
pluto = planet(m=1.48e15, r_body=6.2e3, r_orbit=39.52*AU, phi_solar=0.9)

# Юпитер
jupiter = planet(m=1.8987e27, r_body=69911e3, r_orbit=5.204*AU)
"""

# Локальная гравитационная постоянная
def g(G,M,R):
    return G*M/R**2

# Гравитационный потенциал
def phi_gw(G,M,g,R):
    return G*M/(g*R)

# Вторая космическая скорость (m/s)
def v_esc(G,M,R):
	return math.sqrt(2*G*M/R)

# Расстояние от Земли до космического аппарата(КА)(забыл термин, если знаете -> tarkov3108@gmail.com , буду, фактически, очень признателен)
def l_e(d_s_e,M_s,M_e):
	return d_s_e/(1 + math.sqrt(M_s/M_e))

# Большая полуось
def a(hp, ha):
	return R + (hp + ha)/2

# Орбитальный период(по 3-му закону Кеплера), T = 2pi * sqrt(a^3/mu)
def T(a,G,M):
	return 2*math.pi*math.sqrt(a**3/(G*M))

# Большая полуось (задан T) (m)
def a_T(G,M,T):
	return (G*M*(T/(2*math.pi))**2)**(1./3)

# Орбитальная скорость для эллептического движения
def v(G,M,r,a):
	return math.sqrt(2*G*M/r - G*M/a)

# Орбитальная скорость для гиперболического движения
def v_hyperbolic(G,M,r,a):
	return math.sqrt(2*G*M/r + G*M/a)

# Гравитационный потенциал
def V(G,M,R):
    return G*M/R

# Гомановский переход, deltav1
def deltav1_ht(G,M,r1,r2):
	mu = G*M
	return math.sqrt( (2*mu*r2) / (r1*(r1 + r2)) ) - math.sqrt(mu/r1)

# Гомановский переход, deltav2
def deltav2_ht(G,M,r1,r2):
	mu = G*M
	return - math.sqrt( (2*mu*r1) / (r2*(r1 + r2)) ) + math.sqrt(mu/r2)

# Deltav, необходимый для изменения орбитальной плоскости при пересечении экватора
def deltav_op(Vi,alpha):
	return 2*Vi*math.sin(alpha/2)

# Прецессия линии узлов для низкой околоземной орбиты (в градусах на среднесолнечный день)
def omegap_degperday(i,a,e):
	i = i*180/math.pi
	a = a*1e-3
	return -2.06474e14*math.cos(i)/(a**3.5*(1-e**2)**2)

# Прецессия линии узлов для низкой околоземной орбиты (rad/s) 
def omegap(R,a,e,J2,omega,i):
#	return -3./2. * R**2 / (a * (1. - e**2) )**2 * J2 * omega * math.cos(i)
	return - 3./2 * R**2 / a**2 * J2 * omega * math.cos(i)

# Гиперболическая избыточная скорость (given departure velocity)
def v_inf(G,M,rd,vd):
	return math.sqrt( vd**2 - 2*G*M/rd )

# Вторая космическая(избыточная задана)
def v_d(G,M,rd,v_inf):
	return math.sqrt( 2*G*M/rd + v_inf**2)

# Скорость с гелеоцентрической орбиты (Гравитационный маневр)
def V_Sd(Vp,v_dinf,delta):
	return math.sqrt( Vp**2 + v_dinf**2 - 2*Vp*v_dinf*math.cos(delta) )

# Deltav торможения
def deltav_brake(G,M,vainf,rp,ai):
	return math.sqrt( 2*G*M/rp + vainf**2 ) - math.sqrt( 2*G*M/rp - G*M/ai )

# Уравнение Циалковского
def deltav_rocket(g,Isp,mi,mf):
	return g*Isp*math.log(mi/mf)

# Масса топлива, задана начальной массой космического аппарата (kg)
def m_p(mi,deltav,g0,Isp):
	return mi*(1 - math.exp(-deltav/(g0*Isp)))

# Маневр ориентации по требуемой theta
def t_b_theta(theta_m,n,F,L,I):
	return math.sqrt(theta_m*I/(n*F*L))

# Маневр ориентации по требуемой omega
def t_b_omega(omega_max,n,F,L,I):
	return omega_max*I/(n*F*L)

# Маневр ориентации - использовано топливо
def m_p(n,F,t_b,g,Isp):
	return 2*n*F*t_b/(g*Isp)

# Индуцированное напряжение от троса (закон Фарадея об индукции) при предположении перпендикулярности
def U_tether(V,B,L):
	return V*B*L 

# Вероятность нормальной работы
def R(Lambda,t):
	return math.exp(-Lambda*t)
