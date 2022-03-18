import Module from "./gauss.cjs";

const gauss = await Module();

function scopeCreate() {
    return {
        context: [],
    };
}

function scopeDestroy(scope) {
	  for(let a of scope.context) {
        a.delete();
    }
}

function number(scope, v) {
	let t = gauss.numberFromDouble(v);

	scope.context.push(t);

  return t;
}

function symbol(scope, x) {
	let t = gauss.symbol(x);

	scope.context.push(t);

  return t;
}

function add(scope, a, b) {
    let t = gauss.add(a, b);

    scope.context.push(t);

    return t;
}

let scope = scopeCreate();

let a = number(scope, 4);

let b = number(scope, 0.5);

let c = number(scope, 0.6);

let d = number(scope, 0.333333333333);

let e = symbol(scope, "x");

let f = add(scope, a, b);

let g = add(scope, f, e);

let h = add(scope, g, d);

console.log(gauss.toString(h));

scopeDestroy(scope);

gauss.doLeakCheck();
