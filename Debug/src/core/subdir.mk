################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/core/Angle.o \
../src/core/CS2IO.o \
../src/core/LoopTK.o \
../src/core/PAtom.o \
../src/core/PAtomShell.o \
../src/core/PBlock.o \
../src/core/PBlockConnection.o \
../src/core/PBlockReconnector.o \
../src/core/PBond.o \
../src/core/PCCDSolver.o \
../src/core/PChain.o \
../src/core/PChainNavigator.o \
../src/core/PCluster.o \
../src/core/PConfSpaceNavigator.o \
../src/core/PConformationSpace.o \
../src/core/PDBIO.o \
../src/core/PDataCollector.o \
../src/core/PEnergy.o \
../src/core/PExactIKSolver.o \
../src/core/PGrid.o \
../src/core/PHydrogenBondTracker.o \
../src/core/PInit.o \
../src/core/PMath.o \
../src/core/PMoveAtom.o \
../src/core/PNumRoutines.o \
../src/core/POptimize.o \
../src/core/PPCA.o \
../src/core/PPhiPsiDistribution.o \
../src/core/PProtein.o \
../src/core/PProteinCCDSolver.o \
../src/core/PProteinResidue.o \
../src/core/PResidue.o \
../src/core/PResources.o \
../src/core/PSampMethods.o \
../src/core/PTools.o \
../src/core/PUtilities.o 

CC_SRCS += \
../src/core/Angle.cc \
../src/core/CS2IO.cc \
../src/core/LoopTK.cc \
../src/core/PAtom.cc \
../src/core/PAtomShell.cc \
../src/core/PBlock.cc \
../src/core/PBlockConnection.cc \
../src/core/PBlockReconnector.cc \
../src/core/PBond.cc \
../src/core/PCCDSolver.cc \
../src/core/PChain.cc \
../src/core/PChainNavigator.cc \
../src/core/PCluster.cc \
../src/core/PConfSpaceNavigator.cc \
../src/core/PConformationSpace.cc \
../src/core/PDBIO.cc \
../src/core/PDataCollector.cc \
../src/core/PEnergy.cc \
../src/core/PExactIKSolver.cc \
../src/core/PGrid.cc \
../src/core/PHydrogenBondTracker.cc \
../src/core/PInit.cc \
../src/core/PMath.cc \
../src/core/PMoveAtom.cc \
../src/core/PNumRoutines.cc \
../src/core/POptimize.cc \
../src/core/PPCA.cc \
../src/core/PPhiPsiDistribution.cc \
../src/core/PProtein.cc \
../src/core/PProteinCCDSolver.cc \
../src/core/PProteinResidue.cc \
../src/core/PResidue.cc \
../src/core/PResources.cc \
../src/core/PSampMethods.cc \
../src/core/PTools.cc \
../src/core/PUtilities.cc 

OBJS += \
./src/core/Angle.o \
./src/core/CS2IO.o \
./src/core/LoopTK.o \
./src/core/PAtom.o \
./src/core/PAtomShell.o \
./src/core/PBlock.o \
./src/core/PBlockConnection.o \
./src/core/PBlockReconnector.o \
./src/core/PBond.o \
./src/core/PCCDSolver.o \
./src/core/PChain.o \
./src/core/PChainNavigator.o \
./src/core/PCluster.o \
./src/core/PConfSpaceNavigator.o \
./src/core/PConformationSpace.o \
./src/core/PDBIO.o \
./src/core/PDataCollector.o \
./src/core/PEnergy.o \
./src/core/PExactIKSolver.o \
./src/core/PGrid.o \
./src/core/PHydrogenBondTracker.o \
./src/core/PInit.o \
./src/core/PMath.o \
./src/core/PMoveAtom.o \
./src/core/PNumRoutines.o \
./src/core/POptimize.o \
./src/core/PPCA.o \
./src/core/PPhiPsiDistribution.o \
./src/core/PProtein.o \
./src/core/PProteinCCDSolver.o \
./src/core/PProteinResidue.o \
./src/core/PResidue.o \
./src/core/PResources.o \
./src/core/PSampMethods.o \
./src/core/PTools.o \
./src/core/PUtilities.o 

CC_DEPS += \
./src/core/Angle.d \
./src/core/CS2IO.d \
./src/core/LoopTK.d \
./src/core/PAtom.d \
./src/core/PAtomShell.d \
./src/core/PBlock.d \
./src/core/PBlockConnection.d \
./src/core/PBlockReconnector.d \
./src/core/PBond.d \
./src/core/PCCDSolver.d \
./src/core/PChain.d \
./src/core/PChainNavigator.d \
./src/core/PCluster.d \
./src/core/PConfSpaceNavigator.d \
./src/core/PConformationSpace.d \
./src/core/PDBIO.d \
./src/core/PDataCollector.d \
./src/core/PEnergy.d \
./src/core/PExactIKSolver.d \
./src/core/PGrid.d \
./src/core/PHydrogenBondTracker.d \
./src/core/PInit.d \
./src/core/PMath.d \
./src/core/PMoveAtom.d \
./src/core/PNumRoutines.d \
./src/core/POptimize.d \
./src/core/PPCA.d \
./src/core/PPhiPsiDistribution.d \
./src/core/PProtein.d \
./src/core/PProteinCCDSolver.d \
./src/core/PProteinResidue.d \
./src/core/PResidue.d \
./src/core/PResources.d \
./src/core/PSampMethods.d \
./src/core/PTools.d \
./src/core/PUtilities.d 


# Each subdirectory must supply rules for building sources it contributes
src/core/%.o: ../src/core/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


