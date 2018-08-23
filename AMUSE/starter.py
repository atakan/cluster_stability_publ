import math
import time
import random
import numpy
import sys 
import os

from optparse import OptionParser

from amuse.lab import *
from amuse.datamodel import Particles

from amuse.units import nbody_system , quantities , units, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits

from amuse.io import write_set_to_file

from amuse.ext.halogen_model import new_halogen_model
from amuse.community.halogen.interface import Halogen

output = open('/home/cagri/Desktop/AMUSE/Atakan_Proje/BH_Tree/radiis.txt', 'w')  ## output file
ru_output = open('/home/cagri/Desktop/AMUSE/Atakan_Proje/BH_Tree/r_u_plane.txt', 'w') ## ruplane outputs
rminmax_output = open('/home/cagri/Desktop/AMUSE/Atakan_Proje/BH_Tree/rminmax.txt', 'w') ## spesific star trajectories

spesific_keys = [] # to follow star trajectories

def get_trajectories(stars , clock):
	global rminmax_output
	N = len(stars)
	rminmax_output.write(str(clock.value_in(nbody_system.time)) + "\n")

	for key in spesific_keys:
		j = 0
		while j < N:
			if stars[j].key == key:
				break
			j+=1
		distance = math.sqrt(math.pow(stars[j].x.value_in(nbody_system.length),2) + math.pow(stars[j].y.value_in(nbody_system.length),2) + math.pow(stars[j].z.value_in(nbody_system.length),2) )
		rminmax_output.write(str(distance) + "\n")


def ru_plane_control(stars , clock):
	global ru_output
	N = len(stars)
	i=0
	while i<N:
		distance = math.sqrt(math.pow(stars[i].x.value_in(nbody_system.length),2) + math.pow(stars[i].y.value_in(nbody_system.length),2) + math.pow(stars[i].z.value_in(nbody_system.length),2) )
		radial_velocity = stars[i].vx.value_in(nbody_system.length / nbody_system.time) * stars[i].x.value_in(nbody_system.length) / distance + stars[i].vy.value_in(nbody_system.length / nbody_system.time) * stars[i].y.value_in(nbody_system.length) / distance + stars[i].vz.value_in(nbody_system.length / nbody_system.time) * stars[i].z.value_in(nbody_system.length) / distance
		ru_output.write(str(distance) + "\n")
		ru_output.write(str(radial_velocity) + "\n")
		i+=1
	ru_output.write(str(0) + "\n")

def get_lagrange_radiis(stars , clock):
	global output
	distances = []
	N = len(stars)
	#print N
	i=0
	while i<N:
		distance = math.sqrt(math.pow(stars[i].x.value_in(nbody_system.length),2) + math.pow(stars[i].y.value_in(nbody_system.length),2) + math.pow(stars[i].z.value_in(nbody_system.length),2) )
		distances.append(distance)
		i+=1
	distances = numpy.array(distances)
	distances_sorted = numpy.sort(distances)
	#print distances_sorted
	output.write(str(clock.value_in(nbody_system.time)) + "\n")
	for num in range(1,10):
		output.write(str(distances_sorted[num*N/10]) + "\n")



def initiate_system():

	global spesific_keys

	instance = Halogen(unit_converter=None, redirection='null')
	instance.parameters.convert_nbody = None
	instance.parameters.do_scale = False
	instance.parameters.number_of_particles = 1000
	instance.parameters.alpha = 1
	instance.parameters.beta = 4
	instance.parameters.gamma = 1
	#instance.parameters.black_hole_mass = 1.0 | nbody_system.mass
	instance.parameters.random_seed = random.randint(1, 1e+9)

	#print instance.parameters

	instance.generate_particles()
	result = instance.particles.copy()
	instance.stop()

	########################  --> to obtain spesific stars' keys to follow pathes (for 10 stars)
	i=0
	while i<10:
		spesific_keys.append(result[i].key)
		i+=1
	spesific_keys = numpy.array(spesific_keys)
	########################

	#result.move_to_center()
	return result


def evaluate_system():

	global particles

	while 1:  # Provide VIRIAL THEOREM
		gravity = BHTree() # GRAVITY SOLVER
		gravity.particles.add_particles(particles)
		gravity.parameters.epsilon_squared = 1e-4 | nbody_system.length * nbody_system.length # Softening parameter
		gravity.parameters.opening_angle = 0.5 # 0.01 for black hole existence
		gravity.parameters.timestep = 1e-2 | nbody_system.time # 1e-4 for black hole existence

		if gravity.potential_energy < -1.9 * gravity.kinetic_energy:
			print "Virial provided"
			break
		else:
			gravity.stop()
			particles = initiate_system()
			print "Virial failed"

	channel_from_gravity_to_framework = gravity.particles.new_channel_to(particles) # Channel for retrieving data

	ekin_init = gravity.kinetic_energy
	epot_init = gravity.potential_energy
	print "initial kinetic energy is " , ekin_init
	print "initial potential energy is " , epot_init
	print "initial total energy is " , ekin_init + epot_init

	start_time = time.time()

	end_time = 100
	delta_t = 1e-2
	timerange = nbody_system.time(numpy.arange(0, end_time, delta_t)) # times
	#########################################
	#sys.path.append(os.path.abspath("/home/cagri/Desktop/AMUSE/Atakan_Proje/BH_Tree"))
	#from controller import get_lagrange_radiis

	energy_controller = 0
	for t in timerange:
		gravity.evolve_model(t)
		channel_from_gravity_to_framework.copy()
		get_lagrange_radiis(particles, t)
		get_trajectories(particles, t)

		print gravity.model_time
		#print particles.center_of_mass()
		energy_controller += 1
		if energy_controller == (end_time / delta_t)/20 :  # 20 pieces of cake !
			ru_plane_control(particles, t)
			print "Error in total energy is" ,( (gravity.kinetic_energy + gravity.potential_energy) - (ekin_init + epot_init) ) / (epot_init + ekin_init)
			energy_controller = 0

	#########################################

	print("Elapsed time --> " + str(time.time() - start_time) + "sec.") # Elapsed time of the code

	ekin = gravity.kinetic_energy
	epot = gravity.potential_energy
	print "final total energy is " , ekin + epot

	gravity.stop()
	print "Error is" ,( (ekin + epot) - (ekin_init + epot_init) ) / (epot_init + ekin_init)


#############################################################################################
particles = initiate_system()
#print particles
#print particles.center_of_mass()
#print dynamical_timescale(particles, mass_fraction=0.7)
#particles.scale_to_standard()
#sparticles.move_to_center()

evaluate_system()

#############################################################################################
output.close()
ru_output.close()
rminmax_output.close()