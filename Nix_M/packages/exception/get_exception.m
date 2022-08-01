function [ exception ] = get_exception( type, message )
% 返回NIX:TYPE类型的异常
% type: EMPTY, WARN, ERROR
errID=['Nix_M:' upper(type)];
message=['[' upper(type) '] ' message];
exception=MException(errID,message);
end