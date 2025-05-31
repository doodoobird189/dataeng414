
# make this more generic, e.g.
# find_programs(sox find python)

# find sox
find_program(SOX sox)
if(SOX)
  message(STATUS "Found sox: ${SOX}")
else()
  message(FATAL_ERROR "Can't find sox")
endif()

# find find
find_program(FIND find)
if(FIND)
  message(STATUS "Found find: ${FIND}")
else()
  message(FATAL_ERROR "Can't find find")
endif()

# find python
find_program(PYTHON python)
if(PYTHON)
  message(STATUS "Found python: ${PYTHON}")
else()
  message(FATAL_ERROR "Can't find python")
endif()

function(uniquify_n input output)
  add_custom_command(
    OUTPUT  ${output}
    DEPENDS ${input}
    COMMAND cat `cat ${input}` | sort -u > ${output}
    )
endfunction(uniquify_n)

function(find_target input output path pattern)
  add_custom_command(
    OUTPUT ${output}
    COMMAND ${FIND} ${path} -name "${pattern}" | sort > ${output}
    DEPENDS ${input}
    )
  add_custom_target(
    find-${output} ALL
    DEPENDS ${output}
    )
endfunction(find_target)

