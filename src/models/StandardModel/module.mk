DIR          := models/StandardModel
MODNAME      := StandardModel
SARAH_MODEL  := SM

StandardModel_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

StandardModel_MK     := \
		$(DIR)/module.mk

StandardModel_TWO_SCALE_SUSY_MK := \
		$(DIR)/two_scale_susy.mk

StandardModel_TWO_SCALE_SOFT_MK := \
		$(DIR)/two_scale_soft.mk

StandardModel_TWO_SCALE_MK := \
		$(StandardModel_TWO_SCALE_SUSY_MK) \
		$(StandardModel_TWO_SCALE_SOFT_MK)

StandardModel_SLHA_INPUT := \


StandardModel_GNUPLOT := \
		$(DIR)/StandardModel_plot_rgflow.gnuplot \
		$(DIR)/StandardModel_plot_spectrum.gnuplot

StandardModel_TARBALL := \
		$(MODNAME).tar.gz

LIBStandardModel_SRC :=
EXEStandardModel_SRC :=

LIBStandardModel_HDR :=

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBStandardModel_SRC += \
		$(DIR)/StandardModel_info.cpp \
		$(DIR)/StandardModel_input_parameters.cpp \
		$(DIR)/StandardModel_slha_io.cpp \
		$(DIR)/StandardModel_physical.cpp \
		$(DIR)/StandardModel_utilities.cpp \
		$(DIR)/StandardModel_two_scale_convergence_tester.cpp \
		$(DIR)/StandardModel_two_scale_high_scale_constraint.cpp \
		$(DIR)/StandardModel_two_scale_initial_guesser.cpp \
		$(DIR)/StandardModel_two_scale_low_scale_constraint.cpp \
		$(DIR)/StandardModel_two_scale_model.cpp \
		$(DIR)/StandardModel_two_scale_model_slha.cpp \
		$(DIR)/StandardModel_two_scale_susy_parameters.cpp \
		$(DIR)/StandardModel_two_scale_soft_parameters.cpp \
		$(DIR)/StandardModel_two_scale_susy_scale_constraint.cpp
EXEStandardModel_SRC += \
		$(DIR)/run_StandardModel.cpp \
		$(DIR)/run_cmd_line_StandardModel.cpp \
		$(DIR)/scan_StandardModel.cpp
LIBStandardModel_HDR += \
		$(DIR)/StandardModel_convergence_tester.hpp \
		$(DIR)/StandardModel_high_scale_constraint.hpp \
		$(DIR)/StandardModel_info.hpp \
		$(DIR)/StandardModel_initial_guesser.hpp \
		$(DIR)/StandardModel_input_parameters.hpp \
		$(DIR)/StandardModel_low_scale_constraint.hpp \
		$(DIR)/StandardModel_model.hpp \
		$(DIR)/StandardModel_model_slha.hpp \
		$(DIR)/StandardModel_physical.hpp \
		$(DIR)/StandardModel_slha_io.hpp \
		$(DIR)/StandardModel_spectrum_generator.hpp \
		$(DIR)/StandardModel_susy_scale_constraint.hpp \
		$(DIR)/StandardModel_utilities.hpp \
		$(DIR)/StandardModel_two_scale_convergence_tester.hpp \
		$(DIR)/StandardModel_two_scale_high_scale_constraint.hpp \
		$(DIR)/StandardModel_two_scale_initial_guesser.hpp \
		$(DIR)/StandardModel_two_scale_low_scale_constraint.hpp \
		$(DIR)/StandardModel_two_scale_model.hpp \
		$(DIR)/StandardModel_two_scale_model_slha.hpp \
		$(DIR)/StandardModel_two_scale_soft_parameters.hpp \
		$(DIR)/StandardModel_two_scale_susy_parameters.hpp \
		$(DIR)/StandardModel_two_scale_susy_scale_constraint.hpp

ifneq ($(MAKECMDGOALS),showbuild)
ifneq ($(MAKECMDGOALS),tag)
ifneq ($(MAKECMDGOALS),release)
ifneq ($(MAKECMDGOALS),doc)
-include $(StandardModel_TWO_SCALE_SUSY_MK)
-include $(StandardModel_TWO_SCALE_SOFT_MK)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
ifneq ($(MAKECMDGOALS),pack-$(MODNAME)-src)
ifeq ($(findstring clean-,$(MAKECMDGOALS)),)
ifeq ($(findstring distclean-,$(MAKECMDGOALS)),)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
$(StandardModel_TWO_SCALE_SUSY_MK): run-metacode-$(MODNAME)
		@$(CONVERT_DOS_PATHS) $@
$(StandardModel_TWO_SCALE_SOFT_MK): run-metacode-$(MODNAME)
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
LIBStandardModel_SRC := $(sort $(LIBStandardModel_SRC))
EXEStandardModel_SRC := $(sort $(EXEStandardModel_SRC))

LIBStandardModel_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBStandardModel_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBStandardModel_SRC)))

EXEStandardModel_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEStandardModel_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEStandardModel_SRC)))

LIBStandardModel_DEP := \
		$(LIBStandardModel_OBJ:.o=.d)

EXEStandardModel_DEP := \
		$(EXEStandardModel_OBJ:.o=.d)

LIBStandardModel     := $(DIR)/lib$(MODNAME)$(LIBEXT)

RUN_StandardModel_OBJ := $(DIR)/run_StandardModel.o
RUN_StandardModel_EXE := $(DIR)/run_StandardModel.x

RUN_CMD_LINE_StandardModel_OBJ := $(DIR)/run_cmd_line_StandardModel.o
RUN_CMD_LINE_StandardModel_EXE := $(DIR)/run_cmd_line_StandardModel.x

SCAN_StandardModel_OBJ := $(DIR)/scan_StandardModel.o
SCAN_StandardModel_EXE := $(DIR)/scan_StandardModel.x

METACODE_STAMP_StandardModel := $(DIR)/00_DELETE_ME_TO_RERUN_METACODE

ifeq ($(ENABLE_META),yes)
SARAH_MODEL_FILES_StandardModel := \
		$(shell $(SARAH_DEP_GEN) $(SARAH_MODEL))
endif

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-src \
		distclean-$(MODNAME) run-metacode-$(MODNAME) \
		pack-$(MODNAME)-src

all-$(MODNAME): $(LIBStandardModel)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(StandardModel_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBStandardModel_SRC) $(StandardModel_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBStandardModel_HDR) $(StandardModel_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEStandardModel_SRC) $(StandardModel_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(StandardModel_MK) $(StandardModel_INSTALL_DIR) -m u=rw,g=r,o=r
		install -m u=rw,g=r,o=r $(StandardModel_TWO_SCALE_MK) $(StandardModel_INSTALL_DIR)
ifneq ($(StandardModel_SLHA_INPUT),)
		install -m u=rw,g=r,o=r $(StandardModel_SLHA_INPUT) $(StandardModel_INSTALL_DIR)
endif
		install -m u=rw,g=r,o=r $(StandardModel_GNUPLOT) $(StandardModel_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBStandardModel_DEP)
		-rm -f $(EXEStandardModel_DEP)

clean-$(MODNAME)-obj:
		-rm -f $(LIBStandardModel_OBJ)
		-rm -f $(EXEStandardModel_OBJ)


clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		-rm -f $(LIBStandardModel)
		-rm -f $(RUN_StandardModel_EXE)
		-rm -f $(RUN_CMD_LINE_StandardModel_EXE)
		-rm -f $(SCAN_StandardModel_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

pack-$(MODNAME)-src:
		tar -czf $(StandardModel_TARBALL) \
		$(LIBStandardModel_SRC) $(LIBStandardModel_HDR) \
		$(EXEStandardModel_SRC) \
		$(StandardModel_MK) $(StandardModel_TWO_SCALE_MK) \
		$(StandardModel_SLHA_INPUT) $(StandardModel_GNUPLOT)

$(LIBStandardModel_SRC) $(LIBStandardModel_HDR) $(EXEStandardModel_SRC) \
: run-metacode-$(MODNAME)
		@true

run-metacode-$(MODNAME): $(METACODE_STAMP_StandardModel)
		@true

ifeq ($(ENABLE_META),yes)
$(METACODE_STAMP_StandardModel): $(DIR)/start.m $(DIR)/FlexibleSUSY.m $(META_SRC) $(TEMPLATES) $(SARAH_MODEL_FILES_StandardModel)
		"$(MATH)" -run "Get[\"$<\"]; Quit[]"
		@touch "$(METACODE_STAMP_StandardModel)"
		@echo "Note: to regenerate StandardModel source files," \
		      "please remove the file "
		@echo "\"$(METACODE_STAMP_StandardModel)\" and run make"
		@echo "---------------------------------"
else
$(METACODE_STAMP_StandardModel):
		@true
endif

$(LIBStandardModel_DEP) $(EXEStandardModel_DEP) $(LIBStandardModel_OBJ) $(EXEStandardModel_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBStandardModel_DEP) $(EXEStandardModel_DEP) $(LIBStandardModel_OBJ) $(EXEStandardModel_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIBStandardModel): $(LIBStandardModel_OBJ)
		$(MAKELIB) $@ $^

$(RUN_StandardModel_EXE): $(RUN_StandardModel_OBJ) $(LIBStandardModel) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(RUN_CMD_LINE_StandardModel_EXE): $(RUN_CMD_LINE_StandardModel_OBJ) $(LIBStandardModel) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

$(SCAN_StandardModel_EXE): $(SCAN_StandardModel_OBJ) $(LIBStandardModel) $(LIBFLEXI) $(LIBLEGACY) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(LAPACKLIBS) $(BLASLIBS) $(FLIBS)

ALLDEP += $(LIBStandardModel_DEP) $(EXEStandardModel_DEP)
ALLSRC += $(LIBStandardModel_SRC) $(EXEStandardModel_SRC)
ALLLIB += $(LIBStandardModel)
ALLEXE += $(RUN_StandardModel_EXE) $(RUN_CMD_LINE_StandardModel_EXE) $(SCAN_StandardModel_EXE)
