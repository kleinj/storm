#pragma once

#include <iostream>
#include <sstream>
#include "storm/utility/NumericalTypes.h"
#include "storm/storage/FlexibleSparseMatrix.h"
#include "storm/storage/SCCInfo.h"

namespace storm {
namespace storage {

class HTMLMatrixExport {
public:
    HTMLMatrixExport(const std::string& filename) : out(filename) {
    }

    void setPrintComplexityOption(bool value) {
        complexity = value;
    }

    void setPrintStatsOption(bool value) {
        stats = value;
    }

    void setSCCInfo(boost::optional<const SCCInfo&> sccInfo) {
        this->sccInfo = sccInfo;
    }

    void resetSCCInfo() {
        this->sccInfo = boost::optional<const SCCInfo&>();
    }

    void printHeader(const std::string& title, const std::string& css = "") {
        out << "<!DOCTYPE html>\n";
        out << "<html lang='en'>\n";
        out << "<head>\n";
        out << "<meta charset='utf-8' />\n";
        out << "<title>";
        printText(title);
        out << "</title>\n";
        if (css != "") {
            out << "<style>\n" << css << "</style>\n";
        }
        out << "<script src='https://wwwtcs.inf.tu-dresden.de/~klein/intern/matrix-slides/slides.js'></script>\n";
        out << "<link rel='stylesheet' href='https://wwwtcs.inf.tu-dresden.de/~klein/intern/matrix-slides/slides.css'></link>\n";
        out << "</head>\n";
        out << "<body onload='init()'>\n";
        out << "<h1>";
        printText(title);
        out << "</h1>\n";
    }

    void printFooter() {
        out << "</body>\n";
        out << "</html>" << std::endl;
    }

    void print(const std::string& rawText) {
        out << rawText;
    }

    void printText(const std::string& text) {
        out << text;
    }

    std::string getStandardCSS() {
        std::stringstream css;
        css << " table { border-collapse: collapse; }\n";
        css << " table, th, td { border: 1px solid black; padding: 0.2em;}\n";
        css << " .dia {background-color: lightblue;}\n";
        css << " .sccA {background-color: #a0d0a0;}\n";
        css << " .sccB {background-color: #e0ffe0;}\n";
        css << " .processed {border: 3px solid darkblue;}\n";
        return css.str();
    }

    template <typename T>
    std::string value(const T& t) {
        if (prettyPrint && !complexity && !stats) {
            return prettyPrintHTML(t);
        }
        std::stringstream s;
        s << t;

        std::string content = s.str();

        if (complexity) {
            if (content == "0") {
                return "·";
            }
            std::size_t length = content.length();
            if (prettyPrint) {
                return std::string("<span title='") + s.str() + "'>" + std::to_string(length) + "</span>";
            } else {
                return std::to_string(length);
            }
        } else if (stats) {
            if (content == "0") {
                return ".";
            }
            return std::string("<span title='") + s.str() + "'>" + storm::utility::NumericalTypes::getStats(t) + "</span>";
        }

        if (content == "0") {
            return "·";
        }
        return content;
    }

    template <typename Matrix, typename ValueType>
    struct MatrixPrinter {
        static void printMatrix(HTMLMatrixExport& htmlExport,
                                const Matrix& matrix,
                                boost::optional<std::vector<ValueType> const&> bVector,
                                boost::optional<std::size_t> processed_until_row,
                                boost::optional<std::size_t> highlight) {

            std::ofstream& out = htmlExport.out;

            out << "<table>\n";
            for (std::size_t row = 0; row < matrix.getRowCount(); row++) {
                out << "<tr class='";
                if (htmlExport.sccInfo) {
                    if (htmlExport.sccInfo.get().getSCC(row) % 2 == 0) {
                        out << "sccA";
                    } else {
                        out << "sccB";
                    }
                }
                if ((highlight && row == highlight.get()) ||
                    (processed_until_row && row <= processed_until_row.get())) {
                    out << " processed";
                }
                out << "'>\n";
                for (std::size_t col = 0; col < matrix.getColumnCount(); col++) {
                    out << "<td title='" << row << "," << col << "'";
                    if (col == row) {
                        out << " class='dia'";
                    }
                    out << ">";
                    out << htmlExport.value(matrix(row,col));
                    out << "</td>\n";
                }
                if (bVector) {
                    out << "<td title='" << row << "," << matrix.getColumnCount() << "'>";
                    out << htmlExport.value(bVector.get()[row]);
                    out << "</td>\n";
                }
                out << "</tr>\n";
            }
            out << "</table>" << std::endl;
        }
    };

    // specialization for FlexibleSparseMatrix
    template <typename ValueType>
    struct MatrixPrinter<storm::storage::FlexibleSparseMatrix<ValueType>, ValueType> {
        typedef typename storm::storage::FlexibleSparseMatrix<ValueType> Matrix;

        static void printMatrix(HTMLMatrixExport& htmlExport,
                                const Matrix& matrix,
                                boost::optional<std::vector<ValueType> const&> bVector,
                                boost::optional<std::size_t> processed_until_row,
                                boost::optional<std::size_t> highlight) {
            std::ofstream& out = htmlExport.out;

            out << "<table>\n";
            for (std::size_t row = 0; row < matrix.getRowCount(); row++) {
                out << "<tr class='";
                if (htmlExport.sccInfo) {
                    if (htmlExport.sccInfo.get().getSCC(row) % 2 == 0) {
                        out << "sccA";
                    } else {
                        out << "sccB";
                    }
                }
                if ((highlight && row == highlight.get()) ||
                    (processed_until_row && row <= processed_until_row.get())) {
                    out << " processed";
                }
                out << "'>\n";
                auto it = matrix.row_begin(row);
                auto end = matrix.row_end(row);
                for (std::size_t col = 0; col < matrix.getColumnCount(); col++) {
                    out << "<td title='" << row << "," << col << "'";
                    if (col == row) {
                        out << " class='dia'";
                    }
                    out << ">";
                    if (it == end) {
                        out << htmlExport.value(matrix.getZeroValue());
                    } else if (it->getColumn() == col) {
                        out << htmlExport.value(it->getValue());
                        ++it;
                    } else {
                        assert(it->getColumn() > col);
                        out << htmlExport.value(matrix.getZeroValue());
                    }
                    out << "</td>\n";
                }
                if (bVector) {
                    out << "<td title='" << row << "," << matrix.getColumnCount() << "'>";
                    out << htmlExport.value(bVector.get()[row]);
                    out << "</td>\n";
                }

                out << "</tr>\n";
            }
            out << "</table>" << std::endl;
        }
    };


    // specialization for SparseMatrix
    template <typename ValueType>
    struct MatrixPrinter<storm::storage::SparseMatrix<ValueType>, ValueType> {
        typedef typename storm::storage::SparseMatrix<ValueType> Matrix;

        static void printMatrix(HTMLMatrixExport& htmlExport,
                                const Matrix& matrix,
                                boost::optional<std::vector<ValueType> const&> bVector,
                                boost::optional<std::size_t> processed_until_row,
                                boost::optional<std::size_t> highlight) {
            std::ofstream& out = htmlExport.out;

            ValueType zeroValue(0);

            out << "<table>\n";
            for (std::size_t row = 0; row < matrix.getRowCount(); row++) {
                out << "<tr class='";
                if (htmlExport.sccInfo) {
                    if (htmlExport.sccInfo.get().getSCC(row) % 2 == 0) {
                        out << "sccA";
                    } else {
                        out << "sccB";
                    }
                }
                if ((highlight && row == highlight.get()) ||
                    (processed_until_row && row <= processed_until_row.get())) {
                    out << " processed";
                }
                out << "'>\n";
                auto it = matrix.begin(row);
                auto end = matrix.end(row);
                for (std::size_t col = 0; col < matrix.getColumnCount(); col++) {
                    out << "<td title='" << row << "," << col << "'";
                    if (col == row) {
                        out << " class='dia'";
                    }
                    out << ">";
                    if (it == end) {
                        out << htmlExport.value(zeroValue);
                    } else if (it->getColumn() == col) {
                        out << htmlExport.value(it->getValue());
                        ++it;
                    } else {
                        assert(it->getColumn() > col);
                        out << htmlExport.value(zeroValue);
                    }
                    out << "</td>\n";
                }
                if (bVector) {
                    out << "<td title='" << row << "," << matrix.getColumnCount() << "'>";
                    out << htmlExport.value(bVector.get()[row]);
                    out << "</td>\n";
                }

                out << "</tr>\n";
            }
            out << "</table>" << std::endl;
        }
    };

    template <typename Matrix>
    void printMatrix(const Matrix& matrix,
                     boost::optional<std::size_t> processed_until_row = boost::none,
                     boost::optional<std::size_t> highlight = boost::none) {
        MatrixPrinter<Matrix, typename Matrix::value_type>::printMatrix(*this, matrix, boost::none, processed_until_row, highlight);
    }

    template <typename Matrix, typename ValueType>
    void printMatrix(const Matrix& matrix,
                     std::vector<ValueType> const& bVector,
                     boost::optional<std::size_t> processed_until_row = boost::none,
                     boost::optional<std::size_t> highlight = boost::none) {
        MatrixPrinter<Matrix, ValueType>::printMatrix(*this, matrix, bVector, processed_until_row, highlight);
    }

    template <typename ValueType>
    void printVector(const std::vector<ValueType>& vector) {
        out << "<table>\n";
        for (std::size_t row = 0; row < vector.size(); row++) {
            if (sccInfo) {
                if (sccInfo.get().getSCC(row) % 2 == 0) {
                    out << "<tr class='sccA'>";
                } else {
                    out << "<tr class='sccB'>";
                }
            } else {
                out << "<tr>";
            }
            out << "<td title='" << row << "'>";
            out << value(vector.at(row));
            out << "</td></tr>\n";
        }
        out << "</table>" << std::endl;
    }

    void slideBegin() {
        out << "<div class='slide' id='slide-" << curSlide << "'>\n";
    }

    void slideEnd() {
        out << "</div>" << std::endl;
        curSlide++;
    }

    void flush() {
        out.flush();
    }

    void close() {
        out.close();
    }

    template <typename T>
    inline std::string prettyPrintHTML(const T& t) {
        std::stringstream s;
        s << t;
        return s.str();
    }

    template<typename Pol, bool AS>
    inline std::string prettyPrintHTML(const carl::RationalFunction<Pol, AS>& rhs) {
        std::stringstream s;
        if (rhs.isConstant()) {
            s << rhs.constantPart();
        } else {
            s << "(" << rhs.nominatorAsPolynomial() << ") <big><b>/</b></big> (" << rhs.denominatorAsPolynomial() << ")";
        }
        return s.str();
    }

private:
    std::ofstream out;
    unsigned int curSlide = 0;
    bool prettyPrint = true;
    bool complexity = true;
    bool stats = false;
    boost::optional<const SCCInfo&> sccInfo;
};

}
}
