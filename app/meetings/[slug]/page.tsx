import { notFound } from "next/navigation"
import { getFileBySlug } from "@/lib/server-markdown-loader"
import { MarkdownContent } from "@/components/markdown-content"
import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { CalendarDays, Clock, MapPin, Users, Video, User, FileText, ExternalLink } from "lucide-react"
import { format } from "date-fns"
import { zhCN } from "date-fns/locale"
import Link from "next/link"
import type { Metadata } from "next"

interface MeetingPageProps {
  params: Promise<{
    slug: string
  }>
}

export async function generateStaticParams() {
  // Static generation - return empty array for dynamic rendering
  return []
}

export async function generateMetadata({ params }: MeetingPageProps): Promise<Metadata> {
  const { slug } = await params
  const meetingFile = getFileBySlug('meetings', slug)

  if (!meetingFile) {
    return {
      title: "会议未找到",
    }
  }

  const meeting = meetingFile.metadata
  return {
    title: meeting.title,
    description: meeting.summary,
  }
}

export default async function MeetingPage({ params }: MeetingPageProps) {
  const { slug } = await params
  const meetingFile = getFileBySlug('meetings', slug)

  if (!meetingFile) {
    notFound()
  }

  const meeting = meetingFile.metadata

  const typeIcons = {
    线上: Video,
    线下: MapPin,
    混合: Users,
  }

  const TypeIcon = typeIcons[meeting.eventType as keyof typeof typeIcons]

  const statusColors = {
    即将开始: "bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200",
    进行中: "bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200",
    已结束: "bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200",
  }

  return (
    <div className="min-h-screen bg-background">
      <div className="container mx-auto px-4 py-8 max-w-4xl">
        {/* Header */}
        <div className="mb-8">
          <Link href="/meetings" className="text-sm text-muted-foreground hover:text-primary mb-4 inline-block">
            ← 返回会议列表
          </Link>

          <div className="flex flex-wrap gap-2 mb-4">
            {meeting.tags.map((tag) => (
              <Badge key={tag} variant="secondary" className="text-xs">
                {tag}
              </Badge>
            ))}
            <Badge className={`text-xs ${statusColors[meeting.status as keyof typeof statusColors]}`}>
              {meeting.status}
            </Badge>
          </div>

          <h1 className="text-3xl font-bold mb-4">{meeting.title}</h1>
          <p className="text-lg text-muted-foreground">{meeting.summary}</p>
        </div>

        {/* Meeting Info */}
        <Card className="mb-8">
          <CardHeader>
            <h2 className="text-xl font-semibold">会议信息</h2>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="flex items-center gap-2">
                <CalendarDays className="h-5 w-5 text-muted-foreground" />
                <span>{format(new Date(meeting.date), "yyyy年MM月dd日", { locale: zhCN })}</span>
              </div>

              <div className="flex items-center gap-2">
                <Clock className="h-5 w-5 text-muted-foreground" />
                <span>{meeting.time}</span>
              </div>

              <div className="flex items-center gap-2">
                <TypeIcon className="h-5 w-5 text-muted-foreground" />
                <span>{meeting.eventType}</span>
              </div>

              <div className="flex items-center gap-2">
                <MapPin className="h-5 w-5 text-muted-foreground" />
                <span>{meeting.location}</span>
              </div>
            </div>

            {meeting.speaker && (
              <div className="flex items-center gap-2">
                <User className="h-5 w-5 text-muted-foreground" />
                <span>主讲人: {meeting.speaker}</span>
              </div>
            )}

            {meeting.capacity && (
              <div className="flex items-center gap-2">
                <Users className="h-5 w-5 text-muted-foreground" />
                <span>
                  容量: {meeting.registered || 0}/{meeting.capacity} 人
                </span>
              </div>
            )}

            {meeting.meetingUrl && (
              <div className="pt-4">
                <Button asChild>
                  <Link href={meeting.meetingUrl} target="_blank" rel="noopener noreferrer">
                    <ExternalLink className="h-4 w-4 mr-2" />
                    加入会议
                  </Link>
                </Button>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Materials */}
        {meeting.materials && meeting.materials.length > 0 && (
          <Card className="mb-8">
            <CardHeader>
              <h2 className="text-xl font-semibold">相关材料</h2>
            </CardHeader>
            <CardContent>
              <div className="space-y-2">
                {meeting.materials.map((material, index) => (
                  <div key={index} className="flex items-center gap-2">
                    <FileText className="h-4 w-4 text-muted-foreground" />
                    <Link
                      href={material}
                      className="text-sm hover:text-primary transition-colors"
                      target="_blank"
                      rel="noopener noreferrer"
                    >
                      {material.split('/').pop() || `材料 ${index + 1}`}
                    </Link>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        )}

        {/* Content */}
        <Card>
          <CardContent className="pt-6">
            <MarkdownContent content={meetingFile.htmlContent} />
          </CardContent>
        </Card>
      </div>
    </div>
  )
}