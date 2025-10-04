function(replace_or_add_flag original_flags old_flag new_flag output_var)
  string(REPLACE "${old_flag}" "${new_flag}" modified_flags "${original_flags}")
  if("${modified_flags}" STREQUAL "${original_flags}")
    # If the flag wasn't replaced, it wasn't present, so we add it
    string(STRIP "${original_flags} ${new_flag}" modified_flags)
  endif()
  set(${output_var} "${modified_flags}" PARENT_SCOPE)
endfunction()

function(add_string_if_not_exists variable value)
  # Convert the variable to a list
  string(REPLACE " " ";" var_list "${${variable}}")
  
  # Check if the value is already in the list
  list(FIND var_list "${value}" index)
  
  if(index EQUAL -1)
    # If the value is not in the list, add it
    set(${variable} "${${variable}} ${value}" PARENT_SCOPE)
  endif()
endfunction()
