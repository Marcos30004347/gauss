#include "String.hpp"

string::string(const char *str) {
  size = 0;

  while (str[size] != '\0') {
    size++;
  }

  data = new char[size + 1];

  for (unsigned long i = 0; i < size; i++) {
    data[i] = str[i];
  }

  data[size] = '\0';
}

string::string(const string &str) {
  size = str.size;

  data = new char[size + 1];

  for (unsigned long i = 0; i < size; i++) {
    data[i] = str.data[i];
  }

  data[size] = '\0';
}

string::string(string &&str) {
  size = str.size;

  data = str.data;

  str.data = 0;
}

string &string::operator=(const string &str) {
  if (data)
    delete[] data;

  size = str.size;

  data = new char[size + 1];

  for (unsigned long i = 0; i < size; i++) {
    data[i] = str.data[i];
  }

  data[size] = '\0';

  return *this;
}

string &string::operator=(string &&str) {
  if (data)
    delete[] data;

  size = str.size;

  data = str.data;

  str.data = 0;

  return *this;
}

string::~string() {
  if (data)
    delete[] data;
}

string string::operator+(const string other) {
  string n = "";

  n.size = this->length() + other.length();

  n.data = new char[this->length() + other.length() + 1];

  unsigned long i = 0;

  for (; i < this->length(); i++) {
    n.data[i] = this->data[i];
  }

  for (; i < n.length(); i++) {
    n.data[i] = other.data[i - this->length()];
  }

  n.data[i] = '\0';

  return n;
}

string &string::operator+=(const string other) {
  string n = "";

  n.size = this->length() + other.length();

  n.data = new char[this->length() + other.length() + 1];

  unsigned long i = 0;

  for (; i < this->length(); i++) {
    n.data[i] = this->data[i];
  }

  for (; i < n.length(); i++) {
    n.data[i] = other.data[i - this->length()];
  }

  n.data[i] = '\0';

  *this = n;

  return *this;
}

bool string::operator==(const string other) {
  if (this->length() != other.length())
    return false;

  for (unsigned long i = 0; i < this->length(); i++) {
    if (this->data[i] != other.data[i]) {
      return false;
    }
  }

  return true;
}

void string::reverse() {
    unsigned int i, temp;

    for (i = 0; i < size/2; i++)
    {
        temp = data[i];
        data[i] = data[size - i - 1];
        data[size - i - 1] = temp;
    }
}

int strcmp(const char *X, const char *Y) {
  while (*X) {
    if (*X != *Y) {
      break;
    }

    X++;
    Y++;
  }

  return *(const unsigned char *)X - *(const unsigned char *)Y;
}

string to_string(long long i) {
  string n = "";

  if (i < 0)
    n += "-";

  char k[2] = {'\0', '\0'};

  while (i) {
    k[0] = '0' + (i % 10);

    n += string(k);

    i = i / 10;
  }

	n.reverse();

  return n;
}
string to_string(unsigned long long i) {
  string n = "";

  if (i < 0)
    n += "-";

  char k[2] = {'\0', '\0'};

  while (i) {
    k[0] = '0' + (i % 10);

    n += string(k);

    i = i / 10;
  }

	n.reverse();

  return n;
}

string to_string(int i) {
  string n = "";

  char k[2] = {'\0', '\0'};

  while (i) {
    k[0] = '0' + (i % 10);

    n += string(k);

    i = i / 10;
  }

	n.reverse();

  return n;
}

string to_string(unsigned int i) {
  string n = "";

  char k[2] = {'\0', '\0'};

  while (i) {
    k[0] = '0' + (i % 10);

    n += string(k);

    i = i / 10;
  }

	n.reverse();

  return n;
}

string operator+(const char *a, string b) {
	return string(a) + b;
}
