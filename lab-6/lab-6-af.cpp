
using namespace std;


class Solution {
public:
    int findRotateSteps(string ring, string key) {
        // Получаем длины строк кольца и ключа
        int keyLength = key.size();
        int ringLength = ring.size();

        // Сохраняем позиции каждой буквы кольца (всего 26 букв)
        vector<int> position[26];
        // Асимптотика: O(n), где n = ringLength
        for (int i = 0; i < ringLength; ++i) {
            position[ring[i] - 'a'].push_back(i);
        }

        // Создаем таблицу динамического программирования (dpTable)
        // Размер: keyLength x ringLength
        // Память: O(m * n), где m = keyLength, n = ringLength
        int dpTable[keyLength][ringLength];
        memset(dpTable, 0x3f, sizeof(dpTable)); // Заполняем большим числом (инфиниум)

        // Инициализация первой строки таблицы DP
        // Асимптотика: O(p1), где p1 - количество позиций первой буквы ключа в кольце
        for (int index : position[key[0] - 'a']) {
            dpTable[0][index] = min(index, ringLength - index) + 1;
        }

        // Заполняем DP таблицу
        // Внешний цикл: O(m), где m = keyLength
        for (int i = 1; i < keyLength; ++i) {
            // Перебор всех позиций текущей буквы ключа
            // Асимптотика: O(p2), где p2 - количество позиций текущей буквы ключа в кольце
            for (int j : position[key[i] - 'a']) {
                // Перебор всех позиций предыдущей буквы ключа
                // Асимптотика: O(p3), где p3 - количество позиций предыдущей буквы ключа в кольце
                for (int k : position[key[i - 1] - 'a']) {
                    // Считаем минимальное расстояние для перехода между позициями j и k
                    int stepDiff = min(abs(j - k), ringLength - abs(j - k)) + 1;
                    dpTable[i][j] = min(dpTable[i][j], dpTable[i - 1][k] + stepDiff);
                }
            }
        }

        // Поиск минимального результата для последней буквы ключа
        // Асимптотика: O(p4), где p4 - количество позиций последней буквы ключа в кольце
        int minSteps = INT_MAX;
        for (int index : position[key[keyLength - 1] - 'a']) {
            minSteps = min(minSteps, dpTable[keyLength - 1][index]);
        }

        return minSteps; // Возвращаем минимальное количество шагов
    }
};