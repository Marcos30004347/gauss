// if (typeof fetch === 'undefined') {
//   await import('path').then(path => globalThis.__dirname = path.dirname(import.meta.url));
//   await import('module').then(module => globalThis.require = module.createRequire(import.meta.url));
// }


import Module from "./gauss.cjs";

const gauss = await Module();

console.log(gauss);

function number(scope, v) {
	let t = gauss.numberFromDouble(v);

	scope.push(t);

  return t;
}

function add(scope, a, b) {
    let t = gauss.add(a, b);

    scope.push(t);

    return t;
}

function scopeCreate() {
    return [];
}

function scopeDestroy(scope) {
	  for(let a of scope) {
        a.delete();
    }
}

let scope = scopeCreate();

let a = number(scope, 0.4);

let b = number(scope, 0.5);

let c = number(scope, 0.6);

let d = number(scope, 0.333333333333);

console.log(gauss.toString(a));
console.log(gauss.toString(b));
console.log(gauss.toString(c));
console.log(gauss.toString(d));

let e = add(scope, a, b);

console.log(gauss.toString(e));

scopeDestroy(scope);

gauss.doLeakCheck();
