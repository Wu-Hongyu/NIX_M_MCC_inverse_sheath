function [ exception ] = get_exception( type, message )
% ����NIX:TYPE���͵��쳣
% type: EMPTY, WARN, ERROR
errID=['Nix_M:' upper(type)];
message=['[' upper(type) '] ' message];
exception=MException(errID,message);
end