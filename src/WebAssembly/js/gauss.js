import gaussjs from "gaussjs.wasm.js"

let gauss = {
    sym       = gaussjs.Module.sym,
    num       = gaussjs.Module.num,
    add       = gaussjs.Module.add,
    mul       = gaussjs.Module.mul,
    div       = gaussjs.Module.div,
    pow       = gaussjs.Module.pow,
    fact      = gaussjs.Module.fact,
    inifinity = gaussjs.Module.inifinity,
    reduce    = gaussjs.Module.reduce,
    expand    = gaussjs.Module.expand,
    operand   = gaussjs.Module.operand,
    toString  = gaussjs.Module.toString,
}

export { gauss }