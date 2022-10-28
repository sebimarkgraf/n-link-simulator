import numpy as np
import n_link_sim.double_link
import n_link_sim.quad_link

# We create the array to return the values in here in python for two reasons:
# - We can only return arrays with swig whose size is known before the function call
# - By creating the array in python the interpreter keeps complete control over the object and the data hold by it

def simulate_double_link(states, actions, lengths, masses, inertias, g, friction, dt, dst,
                         use_pd=0, pdSetPoints=np.zeros((4)), pdGain=np.zeros((4))):
    if len(np.shape(states)) == 1:
        states = np.expand_dims(states, 0)
        actions = np.expand_dims(actions, 0)
    num_samples = len(states)
    if dt-1e-3 < 1.0 < dt+1e-3:
        result = np.zeros((num_samples, 2))
    else:
        result = np.zeros((num_samples, 6))
    double_link.simulate(states, actions,  masses, lengths, inertias, dt, g,
                        friction, dst, use_pd, pdSetPoints, pdGain, result)
    return result


def simulate_quad_link(states, actions, lengths, masses, inertias, g, friction, dt, dst):
    if len(np.shape(states)) == 1:
        states = np.expand_dims(states, 0)
        actions = np.expand_dims(actions, 0)
    num_samples = len(states)
    if dt-1e-3 < 1.0 < dt+1e-3:
        result = np.zeros((num_samples,  4))
    else:
        result = np.zeros((num_samples, 12))
    quad_link.simulate(states, actions, masses, lengths, inertias, dt, g, friction, dst, result)
    return result

