#ifndef MATH_STRING
#define MATH_STRING

class string {
private:
  char* data;
	unsigned long size;

public:
	string(const char* str);

	string(const string&);

	string(string&&);

	~string();

	inline unsigned long length() const { return size; }

	string operator+(const string);

	string& operator+=(const string);

	string& operator=(const string&);

	string& operator=(string&&);

	bool operator==(const string other);

	inline char &operator[](unsigned int i) { return data[i]; }

	inline const char* c_str() const { return data; }

	void reverse();

	static string to_string(int i);
	static string to_string(unsigned int i);
	static string to_string(long long i);
	static string to_string(unsigned long long i);
};

int strcmp(const char *X, const char *Y);

string operator+(const char*, string);

#endif
