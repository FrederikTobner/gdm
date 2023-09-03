add_library(${PROJECT_NAME} INTERFACE ${GDM_SOURCES})

# Precompiles dependencies to speed up compilation of the target
if(MSVC)
    # VisualStudio only accepts header files that also have a source file    
    target_precompile_headers(${PROJECT_NAME} INTERFACE pre_compiled_header.hpp pre_compiled_header.cpp) 
else()
    target_precompile_headers(${PROJECT_NAME} INTERFACE pre_compiled_header.hpp)
endif()
