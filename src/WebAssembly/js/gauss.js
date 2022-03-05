import * as Module from "/gaussjs.wasm.js"

console.log(Module)
console.log(Module.toString)
console.log(Module.sym)

let Gauss = Module
// let gauss = {
//     sym       = Module.sym,
//     num       = gaussjs.Module.num,
//     add       = gaussjs.Module.add,
//     mul       = gaussjs.Module.mul,
//     div       = gaussjs.Module.div,
//     pow       = gaussjs.Module.pow,
//     fact      = gaussjs.Module.fact,
//     inifinity = gaussjs.Module.inifinity,
//     reduce    = gaussjs.Module.reduce,
//     expand    = gaussjs.Module.expand,
//     operand   = gaussjs.Module.operand,
//     toString  = gaussjs.Module.toString,
// }

// export { gauss }
export default Gauss;
