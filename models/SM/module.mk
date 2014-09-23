DIR          := models/SM
MODNAME      := SM
SARAH_MODEL  := SM

SM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

SM_MK     := \
		$(DIR)/module.mk

SM_TWO_SCALE_MK := \
		$(DIR)/two_scale_susy.mk \
		$(DIR)/two_scale_soft.mk

SM_SLHA_INPUT := \


SM_GNUPLOT := \
		$(DIR)/SM_plot_rgflow.gnuplot \
		$(DIR)/SM_plot_spectrum.gnuplot

SM_TARBALL := \
		$(MODNAME).tar.gz

LIBSM_SRC :=
EXESM_SRC :=

LIBSM_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBSM_SRC += \
		$(DIR)/SM_info.cpp \
		$(DIR)/SM_slha_io.cpp \
		$(DIR)/SM_physical.cpp \
		$(DIR)/SM_utilities.cpp \
		$(DIR)/SM_two_scale_convergence_tester.cpp \
		$(DIR)/SM_two_scale_high_scale_constraint.cpp \
		$(DIR)/SM_two_scale_initial_guesser.cpp \
		$(DIR)/SM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SM_two_scale_model.cpp \
		$(DIR)/SM_two_scale_susy_parameters.cpp \
		$(DIR)/SM_two_scale_soft_parameters.cpp \
		$(DIR)/SM_two_scale_susy_scale_constraint.cpp
EXESM_SRC += \
		$(DIR)/run_SM.cpp \
		$(DIR)/scan_SM.cpp
LIBSM_HDR += \
		$(DIR)/SM_convergence_tester.hpp \
		$(DIR)/SM_high_scale_constraint.hpp \
		$(DIR)/SM_info.hpp \
		$(DIR)/SM_initial_guesser.hpp \
		$(DIR)/SM_input_parameters.hpp \
		$(DIR)/SM_low_scale_constraint.hpp \
		$(DIR)/SM_model.hpp \
		$(DIR)/SM_physical.hpp \
		$(DIR)/SM_slha_io.hpp \
		$(DIR)/SM_spectrum_generator.hpp \
		$(DIR)/SM_susy_scale_constraint.hpp \
		$(DIR)/SM_utilities.hpp \
		$(DIR)/SM_two_scale_convergence_tester.hpp \
		$(DIR)/SM_two_scale_high_scale_constraint.hpp \
		$(DIR)/SM_two_scale_initial_guesser.hpp \
		$(DIR)/SM_two_scale_low_scale_constraint.hpp \
		$(DIR)/SM_two_scale_model.hpp \
		$(DIR)/SM_two_scale_soft_parameters.hpp \
		$(DIR)/SM_two_scale_susy_parameters.hpp \
		$(DIR)/SM_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(DIR)/two_scale_susy.mk
-include $(DIR)/two_scale_soft.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(DIR)/two_scale_susy.mk: run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(DIR)/two_scale_soft.mk: run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

endif

# remove duplicates in case all algorithms are used
LIBSM_SRC := $(sort $(LIBSM_SRC))
EXESM_SRC := $(sort $(EXESM_SRC))

LIBSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSM_SRC)))

EXESM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXESM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXESM_SRC)))

LIBSM_DEP := \
		$(LIBSM_OBJ:.o=.d)

EXESM_DEP := \
		$(EXESM_OBJ:.o=.d)

LIBSM     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_SM_OBJ := $(DIR)/run_SM.o
RUN_SM_EXE := $(DIR)/run_SM.x

SCAN_SM_OBJ := $(DIR)/scan_SM.o
SCAN_SM_EXE := $(DIR)/scan_SM.x

METACODE_STAMP_SM := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_SM := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBSM)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_SRC) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBSM_HDR) $(SM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXESM_SRC) $(SM_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(SM_MK) $(SM_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(SM_TWO_SCALE_MK) $(SM_INSTALL_DIR)
ifneq ($(SM_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(SM_SLHA_INPUT) $(SM_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(SM_GNUPLOT) $(SM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBSM_DEP)
		-rm -f $(EXESM_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBSM_OBJ)
		-rm -f $(EXESM_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBSM)
		-rm -f $(RUN_SM_EXE)
		-rm -f $(SCAN_SM_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(SM_TARBALL) \
		$(LIBSM_SRC) $(LIBSM_HDR) \
		$(EXESM_SRC) \
		$(SM_MK) $(SM_TWO_SCALE_MK) \
		$(SM_SLHA_INPUT) $(SM_GNUPLOT)

$(LIBSM_SRC) $(LIBSM_HDR) $(EXESM_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_SM)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_SM): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_SM)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_SM)"
		@echo "Note: to regenerate SM source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_SM)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_SM):
		@true
endif

$(LIBSM_DEP) $(EXESM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBSM_DEP) $(EXESM_DEP) $(LIBSM_OBJ) $(EXESM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBSM): $(LIBSM_OBJ)
		$(MAKELIB) $@ $^

$(RUN_SM_EXE): $(RUN_SM_OBJ) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_SM_EXE): $(SCAN_SM_OBJ) $(LIBSM) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBSM_DEP) $(EXESM_DEP)
ALLSRC += $(LIBSM_SRC) $(EXESM_SRC)
ALLLIB += $(LIBSM)
ALLEXE += $(RUN_SM_EXE) $(SCAN_SM_EXE)
