rhea_INC_DIRS := $(shell find $(RHEA_DIR)/include -type d -not -path "*/.svn*")
rhea_INCLUDE  := $(foreach i, $(rhea_INC_DIRS), -I$(i))

libmesh_INCLUDE := $(rhea_INCLUDE) $(libmesh_INCLUDE)

rhea_LIB := $(RHEA_DIR)/librhea-$(METHOD).la

rhea_APP := $(RHEA_DIR)/rhea-$(METHOD)

# source files
rhea_srcfiles    := $(shell find $(RHEA_DIR)/src -name "*.C" -not -name main.C)
rhea_csrcfiles   := $(shell find $(RHEA_DIR)/src -name "*.c")
rhea_fsrcfiles   := $(shell find $(RHEA_DIR)/src -name "*.f")
rhea_f90srcfiles := $(shell find $(RHEA_DIR)/src -name "*.f90")

# object files
rhea_objects := $(patsubst %.C, %.$(obj-suffix), $(rhea_srcfiles))
rhea_objects += $(patsubst %.c, %.$(obj-suffix), $(rhea_csrcfiles))
rhea_objects += $(patsubst %.f, %.$(obj-suffix), $(rhea_fsrcfiles))
rhea_objects += $(patsubst %.f90, %.$(obj-suffix), $(rhea_f90srcfiles))

# plugin files
rhea_plugfiles    := $(shell find $(RHEA_DIR)/plugins/ -name "*.C" 2>/dev/null)
rhea_cplugfiles   := $(shell find $(RHEA_DIR)/plugins/ -name "*.c" 2>/dev/null)
rhea_fplugfiles   := $(shell find $(RHEA_DIR)/plugins/ -name "*.f" 2>/dev/null)
rhea_f90plugfiles := $(shell find $(RHEA_DIR)/plugins/ -name "*.f90" 2>/dev/null)

# plugins
rhea_plugins := $(patsubst %.C, %-$(METHOD).plugin, $(rhea_plugfiles))
rhea_plugins += $(patsubst %.c, %-$(METHOD).plugin, $(rhea_cplugfiles))
rhea_plugins += $(patsubst %.f, %-$(METHOD).plugin, $(rhea_fplugfiles))
rhea_plugins += $(patsubst %.f90, %-$(METHOD).plugin, $(rhea_f90plugfiles))

# rhea main
rhea_main_src    := $(RHEA_DIR)/src/main.C
rhea_app_objects := $(patsubst %.C, %.$(obj-suffix), $(rhea_main_src))

# dependency files
rhea_deps := $(patsubst %.C, %.$(obj-suffix).d, $(rhea_srcfiles)) \
              $(patsubst %.c, %.$(obj-suffix).d, $(rhea_csrcfiles)) \
              $(patsubst %.C, %.$(obj-suffix).d, $(rhea_main_src))

# If building shared libs, make the plugins a dependency, otherwise don't.
ifeq ($(libmesh_shared),yes)
  rhea_plugin_deps := $(rhea_plugins)
else
  rhea_plugin_deps :=
endif

all:: $(rhea_LIB)

$(rhea_LIB): $(rhea_objects) $(rhea_plugin_deps)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(rhea_objects) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(RHEA_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(rhea_LIB) $(RHEA_DIR)

# include RHEA dep files
-include $(rhea_deps)

# how to build RHEA application
ifeq ($(APPLICATION_NAME),rhea)
all:: rhea

rhea: $(rhea_APP)

$(rhea_APP): $(moose_LIB) $(elk_MODULES) $(rhea_LIB) $(rhea_app_objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
          $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(rhea_app_objects) $(rhea_LIB) $(elk_MODULES) $(moose_LIB) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(ADDITIONAL_LIBS)

endif

#
# Maintenance
#
delete_list := $(rhea_APP) $(rhea_LIB) $(RHEA_DIR)/librhea-$(METHOD).*

cleanall:: 
	make -C $(RHEA_DIR) clean 

###############################################################################
# Additional special case targets should be added here
