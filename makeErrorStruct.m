function errorStruct = makeErrorStruct(errorMessage, errorSeverity)
    errorStruct = struct;
    errorStruct.errorMessage = errorMessage;
    errorStruct.errorSeverity = errorSeverity;
end