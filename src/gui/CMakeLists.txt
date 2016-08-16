include(qvtk)

file(GLOB UI_FILES *.ui)
file(GLOB QT_WRAP *.h)
file(GLOB CPP_FILES *.cpp)

qt5_wrap_ui(UISrcs ${UI_FILES})
add_library(gui ${CPP_FILES} ${UISrcs} ${QT_WRAP})
qt5_use_modules(gui Core Gui)
target_link_libraries(gui ${VTK_LIBRARIES} render)
