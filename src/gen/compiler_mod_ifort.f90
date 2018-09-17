! Specifics for ifort
    
module compiler_mod
    use ifport
    implicit none
    ! By including ifport via the compiler_mod, it can be specified in cmake if compiler_mod_ifort or compiler_mod_gfort should be included
    ! ifport is neccessary for ifort to use the system 
        
end module compiler_mod
    