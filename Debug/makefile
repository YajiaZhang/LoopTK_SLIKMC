################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/utils/utils/subdir.mk
-include src/utils/structs/subdir.mk
-include src/utils/math3d/subdir.mk
-include src/utils/math/subdir.mk
-include src/utils/image/subdir.mk
-include src/utils/camera/subdir.mk
-include src/utils/GLdraw/subdir.mk
-include src/utils/subdir.mk
-include src/tests/subdir.mk
-include src/old/subdir.mk
-include src/core/subdir.mk
-include doc/examples/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: LoopTk.exe

# Tool invocations
LoopTk.exe: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cygwin C++ Linker'
	g++  -o "LoopTk.exe" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) LoopTk.exe
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
