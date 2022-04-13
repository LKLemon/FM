#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math

state_pll = {
	"integrator" : 0.0,
	"phaseEst" : 0.0,
	"feedbackI" : 1.0,
	"feedbackQ" : 0.0,
	"ncoOut_0" : 1.0,
	"ncoOut_q_0" : 0.0,
	"trigOffset" : 0
}

def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01, state = state_pll):

	"""
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
	Cp = 2.666
	Ci = 3.555

	# gain for the proportional term
	Kp = (normBandwidth)*Cp
	# gain for the integrator term
	Ki = (normBandwidth*normBandwidth)*Ci

	# output array for the NCO
	ncoOut = np.empty(len(pllIn)+1)
	ncoOut_q = np.empty(len(pllIn)+1)

	# initialize internal state
	integrator = state["integrator"]
	phaseEst = state["phaseEst"]
	feedbackI = state["feedbackI"]
	feedbackQ = state["feedbackQ"]
	ncoOut[0] = state["ncoOut_0"]
	ncoOut_q[0] = state["ncoOut_q_0"]
	trigOffset = state["trigOffset"]
	# note: state saving will be needed for block processing

	for k in range(len(pllIn)):

		# phase detector
		errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

		# four-quadrant arctangent discriminator for phase error detection
		errorD = math.atan2(errorQ, errorI)

		# loop filter
		integrator = integrator + Ki*errorD

		# update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator

		# internal oscillator
		trigOffset += 1
		trigArg = 2*math.pi*(freq/Fs)*(trigOffset) + phaseEst
		feedbackI = math.cos(trigArg)
		feedbackQ = math.sin(trigArg)
		ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
		ncoOut_q[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

	state["integrator"] = integrator
	state["phaseEst"] = phaseEst
	state["feedbackI"] =feedbackI
	state["feedbackQ"] = feedbackQ
	state["ncoOut_0"] = ncoOut[-1]
	state["ncoOut_q_0"] = ncoOut_q[-1]
	state["trigOffset"] = trigOffset

	# for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
	return ncoOut, ncoOut_q, state
	# for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":

	pass
