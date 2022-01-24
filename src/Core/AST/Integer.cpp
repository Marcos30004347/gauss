#include "Integer.hpp"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>


// Int2::Int2() {
// 	this->val = bint<30>::from<int>(0);
// }

// Int2::Int2(const Int2 &a) {
// 	this->val = a.val->copy();
// }

// Int2::Int2(Int2 &&a) {
//   this->val = a.val;
//   a.val = nullptr;
// }

// Int2::Int2(long int v) {
// 	this->val = bint<30>::from<long int>(v);
// }

// Int2::Int2(long long v) {
// 	this->val = bint<30>::from<long long>(v);
// }

// Int2::Int2(unsigned long long v) {
//   this->val = bint<30>::from<unsigned long long>(v);
// }

// Int2::Int2(unsigned int v) {
// 	this->val = bint<30>::from<unsigned int>(v);
// }

// Int2::Int2(int v) {
// 	this->val = bint<30>::from<int>(v);
// }

// Int2::Int2(double b) {
// 	this->val = bint<30>::from(b);
// }

// Int2::~Int2() {
//   if (this->val) delete this->val;
// }

// Int2::Int2(bint<30> *v) {
// 	this->val = v;
// }

// std::string Int2::to_string() { return this->val->to_string(); }


Int::Int() {
  this->flag = 0;
	this->x = 0;
	// assert(this->to_string() == Int2().to_string());
}

Int::Int(const Int &a) {
	if(!a.flag) {
		this->flag = 0;
		this->x = a.x;
		// assert(this->to_string() == Int2(this->x).to_string());
	} else {
		this->flag = 1;
		this->val = a.val->copy();
		// assert(this->to_string() == Int2(this->val->copy()).to_string());
	}
}

Int::Int(Int &&a) {
	if(!a.flag) {
		this->flag = 0;
		this->x = a.x;

		// assert(this->to_string() == Int2(this->x).to_string());
	} else {
		this->flag = 1;
		this->val = a.val;

		// assert(this->to_string() == Int2(this->val->copy()).to_string());
	}

	a.val = nullptr;
}

Int::Int(long int v) {
	if(v < LONG_LONG_MAX) {
		this->flag = 0;
		this->x = v;
	} else {
		this->flag = 1;
		this->val = bint<30>::from<long int>(v);
	}
	// assert(this->to_string() == Int2(v).to_string());
}

Int::Int(long long v) {
	if(v < LONG_LONG_MAX) {
		this->flag = 0;
		this->x = v;
	} else {
		this->flag = 1;
		this->val = bint<30>::from<long long>(v);
	}

	// assert(this->to_string() == Int2(v).to_string());
}

Int::Int(unsigned long long v) {
	if(v < (unsigned long long)LONG_LONG_MAX) {
		this->flag = 0;
		this->x = v;
	} else {
		this->flag = 1;
		this->val = bint<30>::from<unsigned long long>(v);
	}

	// assert(this->to_string() == Int2(v).to_string());
}

Int::Int(unsigned int v) {
	if(v < (unsigned int)LONG_LONG_MAX) {
		this->flag = 0;
		this->x = v;
	} else {
		this->flag = 1;
		this->val = bint<30>::from<unsigned int>(v);
	}

	// assert(this->to_string() == Int2(v).to_string());
}

Int::Int(int v) {
	this->flag = 0;
	this->x = v;

	// assert(this->to_string() == Int2(v).to_string());

}

Int::Int(double b) {
	if(b < (double)LONG_LONG_MAX) {
		this->flag = 0;
		this->x = b;
	} else {
		this->flag = 1;
		this->val = bint<30>::from(b);
	}

	// assert(this->to_string() == Int2(b).to_string());
}

Int::~Int() {
  if (this->flag && this->val) delete this->val;
}

Int::Int(bint<30> *v) {
	this->flag = 1;
	this->val = v;

	// assert(this->to_string() == Int2(v->copy()).to_string());
}

std::string Int::to_string() {
	if(!this->flag) return std::to_string(x);
  return this->val->to_string();
}
