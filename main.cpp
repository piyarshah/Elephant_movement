#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

struct Time {
    string timestampString;
};

struct DataPoint {
    Time timestamp;
    double latitude;
    double longitude;
    double temperature;
};

double custom_stod(const string& str) {
    double result = 0.0;
    int wholePart = 0;
    double fractionalPart = 0.0;
    int divisor = 1;

    bool isNegative = false;
    size_t i = 0;

    if (str[i] == '-') {
        isNegative = true;
        ++i;
    }

    while (i < str.size() && str[i] != '.') {
        wholePart = wholePart * 10 + (str[i] - '0');
        ++i;
    }

    if (i < str.size() && str[i] == '.') {
        ++i;
        while (i < str.size()) {
            fractionalPart = fractionalPart * 10 + (str[i] - '0');
            divisor *= 10;
            ++i;
        }
    }

    result = wholePart + (fractionalPart / divisor);

    if (isNegative) {
        result = -result;
    }

    return result;
}

vector<DataPoint> loadData(const string& filename) {
    vector<DataPoint> data;
    ifstream file(filename);
    string line;

    getline(file, line);

    while (getline(file, line)) {
        DataPoint dp;
        string token;

        size_t pos = line.find(',');
        token = line.substr(0, pos);
        line.erase(0, pos + 1);
        if (token.empty()) {
            continue;
        }
        dp.longitude = custom_stod(token);

        pos = line.find(',');
        token = line.substr(0, pos);
        line.erase(0, pos + 1);
        if (token.empty()) {
            continue;
        }
        dp.latitude = custom_stod(token);

        pos = line.find(',');
        token = line.substr(0, pos);
        line.erase(0, pos + 1);
        if (token.empty()) {
            continue;
        }
        dp.temperature = custom_stod(token);

        dp.timestamp.timestampString = line;

        data.push_back(dp);
    }

    return data;
}

vector<vector<DataPoint>> bucketSortByMonth(const vector<DataPoint>& data) {
    vector<vector<DataPoint>> buckets(12);

    for (const auto& dp : data) {
        string monthStr = dp.timestamp.timestampString.substr(5, 2);
        int month = stoi(monthStr) - 1;
        buckets[month].push_back(dp);
    }

    return buckets;
}

vector<vector<DataPoint>> bucketSortByTemperature(const vector<DataPoint>& data) {
    vector<vector<DataPoint>> buckets(4);

    for (const auto& dp : data) {
        if (dp.temperature >= 20 && dp.temperature < 25) {
            buckets[0].push_back(dp);
        } else if (dp.temperature >= 25 && dp.temperature < 30) {
            buckets[1].push_back(dp);
        } else if (dp.temperature >= 30 && dp.temperature < 35) {
            buckets[2].push_back(dp);
        } else if (dp.temperature >= 35 && dp.temperature < 40) {
            buckets[3].push_back(dp);
        }
    }

    return buckets;
}

void sortDistancesByTime(vector<DataPoint>& bucket) {
    sort(bucket.begin(), bucket.end(), [](const DataPoint& a, const DataPoint& b) {
        string timeA = a.timestamp.timestampString.substr(11);
        string timeB = b.timestamp.timestampString.substr(11);
        return timeA < timeB;
    });
}

double deg2rad(double deg) {
    return (deg * M_PI / 180.0);
}

double haversineDistance(const DataPoint& p1, const DataPoint& p2) {
    const double R = 6371.0;

    double lat1 = deg2rad(p1.latitude);
    double lon1 = deg2rad(p1.longitude);
    double lat2 = deg2rad(p2.latitude);
    double lon2 = deg2rad(p2.longitude);

    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;

    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) *
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

void printAverageDistancesTemperatureWise(const vector<vector<DataPoint>>& temperatureBuckets) {
    static const vector<string> temperatureLabels = {"20-25", "25-30", "30-35", "35-40"};
    for (size_t i = 0; i < temperatureBuckets.size(); ++i) {
        double totalDistance = 0.0;
        size_t totalDataPoints = 0;
        for (const auto& dp : temperatureBuckets[i]) {
            totalDistance += dp.temperature;
            ++totalDataPoints;
        }
        double averageDistance = totalDataPoints == 0 ? 0.0 : totalDistance / totalDataPoints;
        cout << temperatureLabels[i] << " Celsius: Average Distance: " << averageDistance << " km" << endl;
    }
}

void printAverageDistancesAndTemperaturesTimeWise(const vector<vector<double>>& timeBuckets, const vector<vector<double>>& temperatureBuckets) {
    static const vector<string> timeLabels = {"Morning", "Afternoon", "Evening", "Night"};
    for (size_t i = 0; i < timeBuckets.size(); ++i) {
        double totalDistance = 0.0;
        for (size_t j = 0; j < timeBuckets[i].size(); ++j) {
            totalDistance += timeBuckets[i][j];
        }
        double averageDistance = timeBuckets[i].empty() ? 0.0 : totalDistance / timeBuckets[i].size();
        cout << timeLabels[i] << ": Average Distance: " << averageDistance << " km" << endl;

        double totalTemperature = 0.0;
        for (size_t j = 0; j < temperatureBuckets[i].size(); ++j) {
            totalTemperature += temperatureBuckets[i][j];
        }
        double averageTemperature = temperatureBuckets[i].empty() ? 0.0 : totalTemperature / temperatureBuckets[i].size();
        cout << timeLabels[i] << ": Average Temperature: " << averageTemperature << " Celsius" << endl;
    }
}

void calculateAverageDistancesAndTemperaturesTimeWise(const vector<DataPoint>& data) {
    if (data.empty()) {
        cout << "No data available." << endl;
        return;
    }

    string currentMonth = data[0].timestamp.timestampString.substr(5, 2);
    cout << "Month: " << currentMonth << endl;

    vector<vector<double>> timeBuckets(4);
    vector<vector<double>> temperatureBuckets(4);

    for (size_t i = 0; i < data.size() - 1; ++i) {
        string month = data[i].timestamp.timestampString.substr(5, 2);
        if (month != currentMonth) {
            printAverageDistancesAndTemperaturesTimeWise(timeBuckets, temperatureBuckets);
            currentMonth = month;
            cout << "Month: " << currentMonth << endl;
            timeBuckets.assign(4, vector<double>());
            temperatureBuckets.assign(4, vector<double>());
        }
        double distance = haversineDistance(data[i], data[i + 1]);
        int hour = stoi(data[i].timestamp.timestampString.substr(11, 2));
        int bucketIndex = 0;
        if (hour >= 6 && hour < 12) {
            bucketIndex = 0;
        } else if (hour >= 12 && hour < 18) {
            bucketIndex = 1;
        } else if (hour >= 18 && hour < 24) {
            bucketIndex = 2;
        } else {
            bucketIndex = 3;
        }
        timeBuckets[bucketIndex].push_back(distance);
        temperatureBuckets[bucketIndex].push_back(data[i].temperature);
    }
    printAverageDistancesAndTemperaturesTimeWise(timeBuckets, temperatureBuckets);
}

int main() {
    string filename = "/Users/piyashah/Downloads/Elephants.csv";
    vector<DataPoint> data = loadData(filename);

    vector<vector<DataPoint>> buckets = bucketSortByMonth(data);

    for (auto& bucket : buckets) {
        sortDistancesByTime(bucket);
    }

    for (size_t month = 0; month < buckets.size(); ++month) {
        if (!buckets[month].empty()) {
            calculateAverageDistancesAndTemperaturesTimeWise(buckets[month]);
        }
    }

    vector<vector<DataPoint>> temperatureBuckets = bucketSortByTemperature(data);

    cout << "Average Distances for Each Temperature Range:" << endl;
    printAverageDistancesTemperatureWise(temperatureBuckets);

    return 0;
}

