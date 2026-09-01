#include "EDMD.h"

#include <stdbool.h>
#include <stdio.h>

enum { SIMULATION_FAULT_MESSAGE_CAPACITY = 512 };

static bool simulationFaultActive = false;
static bool simulationFaultRecoverableOverlap = false;
static char simulationFaultDescription[SIMULATION_FAULT_MESSAGE_CAPACITY] = {0};

void simulationClearFault(void){
	simulationFaultActive = false;
	simulationFaultRecoverableOverlap = false;
	simulationFaultDescription[0] = '\0';
}

bool simulationHasFault(void){
	return simulationFaultActive;
}

const char* simulationFaultMessage(void){
	return simulationFaultDescription;
}

void simulationReportFault(const char* message){
	if (simulationFaultActive)
		return;
	if (message == NULL || message[0] == '\0')
		message = "The simulation entered an invalid state and must be restarted.";
	(void)snprintf(simulationFaultDescription,
		sizeof(simulationFaultDescription), "%s", message);
	simulationFaultDescription[sizeof(simulationFaultDescription) - 1] = '\0';
	simulationFaultRecoverableOverlap = false;
	simulationFaultActive = true;
}

void simulationReportRecoverableOverlapFault(const char* message){
	if (simulationFaultActive)
		return;
	simulationReportFault(message);
	simulationFaultRecoverableOverlap = simulationFaultActive;
}

bool simulationFaultIsRecoverableOverlap(void){
	return simulationFaultActive && simulationFaultRecoverableOverlap;
}
