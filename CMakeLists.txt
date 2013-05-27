################################################################################
# em_plugin
#
# A UG4 plugin for discretizations of the Maxwell equations.
################################################################################

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	Electromagnetism)
set(pluginPath	plugins/experimental/electromagnetism)
set(SOURCES		em_plugin.cpp)


################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)

# End of File
