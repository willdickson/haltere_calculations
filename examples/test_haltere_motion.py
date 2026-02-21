import functools
import sympy 
from haltere_calculations import haltere_position
from haltere_calculations import haltere_velocity
from haltere_calculations import coriolis_force
from haltere_calculations import normal_coriolis_force

sympy.init_printing()

class CosFunc(sympy.Function):
    @classmethod
    def eval(cls, amp, freq, t):
        return amp*sympy.cos(2*sympy.pi*freq*t)

explicit = True 
verbose = False 

param = {
        't'     : sympy.Symbol('t'),        # time
        'f'     : sympy.Symbol('f'),        # frequency
        'Phi'   : sympy.Symbol('Phi'),      # haltere stroke amplitude
        'w_1'   : sympy.Symbol('omega_1'),  # x-component fly's angular velocity 
        'w_2'   : sympy.Symbol('omega_2'),  # y-component fly's angular velocity
        'w_3'   : sympy.Symbol('omega_3'),  # z-component fly's angular velocity
        'm'     : sympy.Symbol('m'),        # haltere end knob mass
        'L'     : sympy.Symbol('L'),        # haltere length (base to end knob CM)
        'b'     : sympy.Symbol('b'),        # haltere separation (widthwise)
        'beta'  : sympy.Symbol('beta'),     # haltere tilt angle
        }

if explicit:
    theta = functools.partial(CosFunc, param['Phi'], param['f'])
else:
    theta = sympy.Function('theta')
param['theta'] = theta


haltere_position(param, verbose=verbose)
haltere_velocity(param, verbose=verbose)
coriolis_force(param, verbose=verbose)
force_r, force_l, w_contrib = normal_coriolis_force(param, verbose=verbose)


f_pitch = w_contrib['right'][param['w_1']]
print()
print('f pitch = ')
sympy.pprint(f_pitch)
print('-'*50)

f_roll = w_contrib['right'][param['w_2']]
print()
print('f roll')
sympy.pprint(f_roll)
print('-'*50)

f_yaw = w_contrib['right'][param['w_3']]
print()
print('f yaw')
sympy.pprint(f_yaw)
print('-'*50)









