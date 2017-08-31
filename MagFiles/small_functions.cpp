using namespace std;                    // makes std::cout , std::string unnecessary

string double_to_string(double number);
string int_to_string(int number);

string double_to_string(double number) {
	stringstream sstring;
    sstring << number;
	string s;
	sstring >> s;
	return s;
}

string int_to_string(int number) {
	stringstream sstring;
    sstring << number;
	string s;
	sstring >> s;
	return s;
}
