#include "Definitions.h"
#include "Utils.h"

#define CONCAT_STR(str1, str2) (str1 "" str2)




NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);


bool operator==(const double& v, const PointDouble3D& pp)
{
    return (pp == v);
}

bool operator!=(const double &v, const PointDouble3D &pp)
{
    return (pp != v);
}

bool operator==(const double& v, const PointDouble3DLight& pp)
{
    return (pp == v);
}

bool operator!=(const double &v, const PointDouble3DLight &pp)
{
    return (pp != v);
}

PointDouble3D& PointDouble3D::operator=(const PointDouble3DLight& point)
{
    this->x = point.x;
    this->y = point.y;
    this->z = static_cast<double>(point.z);
    return *this;
}


DataContent::DataContent(const char *name)
{
    (void)name;
    content_data = nullptr;
    content_size = 0;
}

DataContent::DataContent(const DataContent& dc)
{
    set_data(dc.data(), dc.size());
}

DataContent::DataContent(DataContent &&_other) noexcept: content_size(_other.content_size), content_data(_other.content_data)
{
    _other.content_data = nullptr;
    _other.content_size = 0;
}

DataContent::~DataContent()
{
    clear();
}

void DataContent::reserve(const size_t &res_size)
{
    clear();
    content_data = new uint8_t[res_size];
    content_size = res_size;
}

void DataContent::set_data(const void* input_data, const size_t &input_size)
{
    clear();
    content_data = new uint8_t[input_size];
    content_size = input_size;
    memcpy(content_data, input_data, content_size);
}

void DataContent::insert_to_beginning(const void* input_data, const size_t &input_size)
{
    size_t new_size = content_size + input_size;

    uint8_t* tmp_buffer = new uint8_t[new_size];
    memcpy(tmp_buffer, input_data, input_size);
    memcpy(tmp_buffer + input_size, content_data, content_size);
    delete [] content_data;
    content_data = tmp_buffer;

    content_size = new_size;
}

void DataContent::add_to_end(const void* input_data, const size_t &input_size)
{
    size_t new_size = content_size + input_size;

    uint8_t* tmp_buffer = new uint8_t[new_size];
    memcpy(tmp_buffer, content_data, content_size);
    memcpy(tmp_buffer + content_size, input_data, input_size);
    delete [] content_data;
    content_data = tmp_buffer;
    content_size = new_size;
}

void DataContent::clear()
{
    if(content_data && content_size > 0)
    {
        delete [] content_data;
        content_data = nullptr;
    }
    content_size = 0;
}
