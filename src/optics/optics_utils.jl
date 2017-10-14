function get_nonlinearcoeff(res, mode, source)
    getω(source)/c * res.material.nonlinearindex * mode.corefraction(source) / mode.effectivearea(source)
end
